#include <iostream>
#include <stdexcept>
#include <string>
//#include <stack>
#include <vector>
#include <functional>

// utility

// Kada pamtimo operatore, pamtimo i gdje je koji operator nadjen da bi kasnije mogla provjeriti da li je prvi
// operand u zagradi npr
struct Operator {
    char op;
    int index;

    Operator(char c, int i) : op(c), index(i) {}
};

typedef std::vector<Operator> Operators;
typedef std::function<std::string(std::string, std::string, bool)> SolverFun;

int abs(int n) {
    if ( n < 0 ) return -n;
    return n;
}

std::string toString(int n) {
    std::string res;
    if ( n == 0 ) return "0";

    bool negative = false;
    if ( n < 0 ) {
        n = -n;
        negative = true;
    }

    while ( n ) {
        res += (n % 10) + '0';
        n /= 10;
    }

    std::string str;
    if ( negative ) {
        str += '-';
    }

    while ( res.size() ) {
        str += res[res.length() - 1];
        res.erase(res.begin() + res.length() - 1);
    }

    return str;
}

bool isNum(char c) {
    return c >= '0' && c <= '9';
}

bool isChar(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

bool isPara(char c) {
    return c == '(' || c == ')';
}

bool isOperator(char c) {
    return (c == '+') ||
           (c == '-') ||
           (c == '*') ||
           (c == '^');
}

bool isSeparator(char c) {
    //faktore razdvaja sve sto nije broj ili slovo
    return !(isNum(c) || isChar(c));
}

void trim(std::string& s) {
    while ( s.length() && s[0] == ' ' ) {
        s.erase(s.begin());
    }
    while ( s.length() && s[s.length() - 1] == ' ' ) {
        s.erase(s.length() - 1);
    }
}

void removeVar(std::string& expr) {
    expr.erase(expr.begin() + expr.size() - 1);
    //ako je bilo -x, otkinuli x, treba ostati -1
    if ( expr == "-" ) {
        expr = "-1";
    } else if ( expr == "" ) {
        expr = "1";
    }
}

int toNumber(std::string str) {
    trim(str);

    int res = 0;
    bool negative;

    if ( str[0] == '-' ) {
        negative = true;
        str.erase(str.begin());
    } else {
        negative = false;
    }

    for ( char c : str ) {
        if ( !isNum(c) ) {
            throw std::invalid_argument("Nije broj!");
        }

        res = (res * 10) + (c - '0');
    }

    if ( negative ) return -res;
    return res;
}

void splitExpression(std::string exp, std::vector<std::string>& operands, Operators& operators);


bool isVariable(std::string var) {
    //token je varijabla ako sadrzi karakter
    trim(var);

    //smije biti negativan
    if ( var[0] == '-' ) {
        var.erase(var.begin());
    }

    for ( char c : var ) {
        if ( isChar(c) ) return true;
    }
    return false;
}

bool isMinimized(const std::string& expr) {
    //maksimalno je sredjen ako nema zagrada i svaki stepen uz x ima maksimalno jednom
    //takodjer ima maksimalno jedan slobodni koeficijent

    if ( expr.find('(') != std::string::npos || expr.find(')') != std::string::npos ) {
        return false;
    }

    std::vector<std::string> operands;
    Operators operators;

    splitExpression(expr, operands, operators);

    std::vector<std::string> exponents;

    for ( auto op : operands ) {
        if ( !isVariable(op)) {
            exponents.push_back("0");
        }
    }
    if ( exponents.size() > 1 ) {
        //imamo vise od jednog slobodnog koeficijenta
        return false;
    }

    for ( int i = 0; i < operators.size(); ++ i ) {
        if ( operators[i].op == '^' ) {
            // stepenovanje se smije javiti akko je uz prvi operand x (uz drugi ne smije biti)
            if ( isVariable(operands[i+1])) {
                throw std::invalid_argument("Exponent ne smije sadrzavati nepoznatu!");
            }
            if ( !isVariable(operands[i]) ) {
                return false; //imamo const ^ const, moze se to jos malo srediti...
            } else {
                //jeste varijabla, pa hajd zapamti kad si vec tu na koji stepe ide.
                for ( auto n : exponents ) {
                    if ( n == operands[i+1] ) {
                        return false; // imamo 2 puta x^exp sa istim exp, pa se moze jos srediti...
                    }
                }
            }
        }
    }

    return true;
}


//vraca broj operanada u nekom izrazu, npr
// 2x + 1 - 6 -> 3
// 2x + x^2 -> 2
// 6x^1 -> 1
int operandsCount(const std::string& str) {
    int res = 1;
    for ( char c : str ) {
        if ( isOperator(c) ) {
            ++res;
        }
    }
    return res;
}

std::string solve_exp(std::string base, std::string exp, bool baseInsidePara);
std::string solve_minus(std::string a, std::string b, bool aInsidePara);
std::string solve_plus(std::string a, std::string b, bool aInsidePara);
std::string solve_mul(std::string a, std::string b, bool aInsidePara);




struct OperationSolver {
    char op;
    SolverFun solver;
};

class Polinom {
private:
    struct Koeficijent {
        int element;
        int exponent;

        Koeficijent(int el, int ex) : element(el), exponent(ex) {}
        Koeficijent(const Koeficijent& k) : element(k.element), exponent(k.exponent) {}
    };

    struct Cvor {
        Koeficijent koeficijent;
        Cvor* next;
        Cvor(int element, int exponent, Cvor* next = nullptr) :
            koeficijent(element, exponent), next(next) {}

        explicit Cvor(const Koeficijent& k, Cvor* next = nullptr) : koeficijent(k), next(next) {}
    };

    Cvor* m_first;

    void copyFrom(const Polinom& src);
    void takeOverFrom(Polinom& src);
    void release();

    Cvor* findByExp(int exp) const;
public:
    Polinom() : m_first(nullptr) {}
    Polinom(const Polinom& src) : m_first(nullptr) {
        copyFrom(src);
    }
    Polinom(Polinom&& src) noexcept : m_first(nullptr) {
        takeOverFrom(src);
    }

    static Polinom ParsirajPolinom(const std::string& exp) {
        // std::cout<<"Parsing -> "<<exp<<std::endl;

        //prvo trivijalne provjere, da li imamo uzastopnih operatora (npr 5--1)
        //ako da, baci izuzetak
        bool found = false;
        for ( char c : exp ) {
            if ( isOperator(c) ) {
                if ( found ) {
                    throw std::invalid_argument("Ima vise uzastopnih operatora!");
                } else {
                    found = true;
                }
            } else {
                found = false;
            }
        }

        //dalje, da li su sve zagrade ispravno otvorene/zatvorene?
        int paraCnt = 0;
        for ( char c : exp ) {
            if ( isPara(c) ) {
                if ( c == ')' ) {
                    --paraCnt;
                } else if ( c == '(' ) {
                    ++paraCnt;
                }

                if ( paraCnt < 0 ) {
                    throw std::invalid_argument("Neispravne zagrade!");
                }
            }
        }
        if ( paraCnt != 0 ) {
            throw std::invalid_argument("Neispravne zagrade!");
        }




        Polinom res;
        Operators operators;
        std::vector<std::string> operands;
        splitExpression(exp, operands, operators);

        if ( !isMinimized(exp) ) {
            //izraz se moze srediti, prvo obavi to pa onda parsiraj

            std::vector<std::vector<OperationSolver>> solve_queue = {
                    {{'^', solve_exp}},
                    {{'*', solve_mul}},
                    {{'+', solve_plus}, {'-', solve_minus}}
            };

            for (const auto &currQueue : solve_queue) {
                for (unsigned int i = 0; i < operators.size(); ++i) {
                    bool exists = false;
                    OperationSolver *solver;

                    for (auto op : currQueue) {
                        if (op.op == operators[i].op) {
                            exists = true;
                            solver = &op;
                            break;
                        }
                    }
                    if (!exists) continue;

                    auto op = operators[i];
                    auto op1 = operands[i];
                    auto op2 = operands[i + 1];
                    //provjeri da li je prije operanda zagrada,
                    // odnosno da li je prvi operand u zagradi
                    bool op1InPara = false;
                    for (int k = op.index; k >= 0; -- k ) {
                        if ( exp[k] == ')' ) {
                            op1InPara = true;
                            break;
                        } else if ( exp[k] == ' ' || exp[k] == op.op ) {
                            continue;
                        } else break;
                    }
                    auto result = solver->solver(op1, op2, op1InPara);

                    // std::cout << "(" << op1 << ")" << op.op << "(" << op2 << ") = " << result << std::endl;

                    operands.erase(operands.begin() + i + 1);
                    operators.erase(operators.begin() + i);

                    operands[i] = result;
                    --i;
                }
            }

        }


        //izraz je sredjen, rastavi koeficijente i sastavi listu
        for ( int i = 0; i < operands.size(); ++ i ) {
            auto& op = operands[i];

            if ( isVariable(op) ) {
                //ubacujem ono sto ide uz varijablu

                //ako je samo jedan clan, npr "2x", prvo skontaj stepen...
                if ( operandsCount(op) == 1 ) {
                    int exponent = 1;
                    if (operators.size() && operators[i].op == '^') {
                        ++i;
                        exponent = toNumber(operands[i]);
                    }
                    //pa onda otkini varijablu...
                    removeVar(op);

                    res.dodajKoeficijent(toNumber(op), exponent);
                } else {
                    //a isto tako moze biti i izraz (npr 1 + 2x^2)
                    //moramo rastaviti na proste clanove, sto vec imamo u funkciji splitExpression
                    int exponent = 1;
                    Operators exprOperators;
                    std::vector<std::string> exprOperands;
                    splitExpression(op, exprOperands, exprOperators);

                    for ( int j = 0; j < exprOperands.size(); ++ j ) {
                        auto& exprOp = exprOperands[j];
                        auto negative = false;

                        if ( exprOperators.size() > j && j > 0 ) {
                            negative = exprOperators[j-1].op == '-';
                        }

                        exponent = 0;

                        if ( isVariable( exprOp ) ) {
                            removeVar(exprOp);

                            if ( exprOperators.size() > j && exprOperators[j].op == '^' ) {
                                exponent = toNumber(exprOperands[j+1]);
                                ++j;
                            } else {
                                exponent = 1;
                            }
                        }
                        if ( negative ) {
                            exprOp = "-" + exprOp;
                        }

                        res.dodajKoeficijent(toNumber(exprOp), exponent);
                    }
                }

            } else {
                //ubacujem slobodni koeficijent
                res.dodajKoeficijent(toNumber(op), 0);
            }
        }
        return res;

    }


    ~Polinom() {
        release();
    }


    Polinom& operator=(const Polinom& src);
    Polinom& operator=(Polinom&& src) noexcept;

    std::string convertToString() const {
        bool first = true;
        std::string result;
        for ( auto it = m_first; it; it = it->next ) {
            if ( it->koeficijent.element == 0 ) {
                continue;
            }

            if ( first ) {
                first = false;
                // da li je prvi clan negativan?
                if ( it->koeficijent.element < 0 ) {
                    result += "-";
                }
            } else {
                if ( it->koeficijent.element < 0 ) {
                    result += " - ";
                } else {
                    result += " + ";
                }
            }

            result += toString(abs(it->koeficijent.element));
            if (it->koeficijent.exponent != 0) {
                result += "x";
                if (it->koeficijent.exponent > 1) {
                    result += "^" + toString(it->koeficijent.exponent);
                }
            }
        }
        return result;
    }
    void Ispisi(std::ostream& stream = std::cout) const;

    Koeficijent* dajKoeficijent(int exponent);
    const Koeficijent* dajKoeficijent(int exponent) const;

    Koeficijent* dodajKoeficijent(int element, int exponent);


    //operatori
    Polinom& operator+=(const Polinom& other);
    Polinom& operator-=(const Polinom& other);
    Polinom& operator*=(const Polinom& other);
};

std::ostream& operator<<(std::ostream& str, const Polinom& p);
std::istream& operator>>(std::istream& str, Polinom& p);

Polinom operator+(const Polinom& a, const Polinom& b);
Polinom operator-(const Polinom& a, const Polinom& b);
Polinom operator*(const Polinom& a, const Polinom& b);


// - Implementacije ispod -

Polinom& Polinom::operator*=(const Polinom& other) {
    //mnozim svaki sa svakim...
    Polinom tmp;
    for ( auto it = m_first; it; it = it->next ) {
        for ( auto it2 = other.m_first; it2; it2 = it2->next ) {
            int exp = it->koeficijent.exponent + it2->koeficijent.exponent;
            int elm = it->koeficijent.element * it2->koeficijent.element;

            if ( auto koef = tmp.dajKoeficijent(exp) ) {
                //imam vec taj exponent, dodaj ga na trenutni
                koef->element += elm;
            } else {
                //novi exponent, ubaci ga na trenutni
                tmp.dodajKoeficijent(elm, exp);
            }
        }
    }

    release();
    takeOverFrom(tmp);

    return *this;
}

Polinom& Polinom::operator+=(const Polinom& other) {
    //idemo kroz drugu listu i za svaki koeficijent provjerimo da li ima u nasoj listi
    //ako ima, inkrementiramo koeficijent, ako ne, dodamo ga.

    for ( auto it = other.m_first; it; it = it->next ) {
        if ( dajKoeficijent(it->koeficijent.exponent) == nullptr ) {
            //nema ovaj stepen, ubaci ga
            dodajKoeficijent(it->koeficijent.element, it->koeficijent.exponent);
        } else {
            //ima ovaj stepen, povecaj ga
            dajKoeficijent(it->koeficijent.exponent)->element += it->koeficijent.element;
        }
    }

    return *this;
}

Polinom& Polinom::operator-=(const Polinom& other) {
    //isto kao += samo sa minus

    for ( auto it = other.m_first; it; it = it->next ) {
        if ( dajKoeficijent(it->koeficijent.exponent) == nullptr ) {
            //nema ovaj stepen, ubaci ga
            dodajKoeficijent(-it->koeficijent.element, it->koeficijent.exponent);
        } else {
            //ima ovaj stepen, smanji ga
            dajKoeficijent(it->koeficijent.exponent)->element -= it->koeficijent.element;
        }
    }

    return *this;
}

Polinom operator+(const Polinom& a, const Polinom& b) {
    Polinom res(a);
    res += b;
    return res;
}

Polinom operator-(const Polinom& a, const Polinom& b) {
    Polinom res(a);
    res -= b;
    return res;
}

Polinom operator*(const Polinom& a, const Polinom& b) {
    Polinom res(a);
    res *= b;
    return res;
}

std::ostream& operator<<(std::ostream& str, const Polinom& p) {
    p.Ispisi(str);
    return str;
}

std::istream &operator>>(std::istream &str, Polinom &p) {
    std::string expression;

    std::getline(str, expression);

    p = Polinom::ParsirajPolinom(expression);
    return str;
}


Polinom &Polinom::operator=(Polinom &&src) noexcept {
    if ( this == &src ) return *this;
    release();
    takeOverFrom(src);
    return *this;
}

Polinom &Polinom::operator=(const Polinom &src) {
    if ( this == &src ) return *this;
    release();
    copyFrom(src);
    return *this;
}

void Polinom::Ispisi(std::ostream& stream) const {
    stream<<convertToString();
}

void Polinom::release() {
    while ( m_first ) {
        auto t = m_first->next;
        delete m_first;
        m_first = t;
    }
}

void Polinom::copyFrom(const Polinom &src) {
    m_first = nullptr;
    Cvor* start = nullptr;

    for ( auto it = src.m_first; it; it = it->next ) {
        if ( !m_first ) {
            m_first = new Cvor(it->koeficijent);
            start = m_first;
        } else {
            m_first->next = new Cvor(it->koeficijent);
            m_first = m_first->next;
        }
    }

    m_first = start;
}

void Polinom::takeOverFrom(Polinom &src) {
    m_first = src.m_first;
    src.m_first = nullptr;
}


Polinom::Koeficijent *Polinom::dajKoeficijent(int exponent) {
    auto res = findByExp(exponent);
    if ( !res ) return nullptr;
    return &res->koeficijent;
}

const Polinom::Koeficijent *Polinom::dajKoeficijent(int exponent) const {
    auto res = findByExp(exponent);
    if ( !res ) return nullptr;
    return &res->koeficijent;
}

Polinom::Cvor *Polinom::findByExp(int exp) const {
    for ( auto it = m_first; it; it = it->next ) {
        if ( it->koeficijent.exponent == exp ) return it;
    }
    return nullptr;
}

Polinom::Koeficijent *Polinom::dodajKoeficijent(int element, int exponent) {
    Cvor* it = m_first;
    Cvor* prev = nullptr;
    while ( it ) {
        if ( it->koeficijent.exponent == exponent ) {
            //ima vec ovaj exponent! problem...
            throw std::invalid_argument("Imam vec taj exponent!");
        } else if ( it ->koeficijent.exponent > exponent ) {
            //nasla sam prvi veci od ovoga, znaci treba ubaciti
            //novi izmedju trenutnog i prethodnog (ako ga ima)
            auto nCvor = new Cvor(element, exponent, it);

            if ( prev ) {
                prev->next = nCvor;
            } else {
                //novi postaje prvi cvor u listi
                m_first = nCvor;
            }

            return &nCvor->koeficijent;
        }

        prev = it;
        it = it->next;
    }

    //ako sam ovdje dosla, treba dodati novi cvor na kraj liste
    auto nCvor = new Cvor(element, exponent);
    if ( !prev ) {
        m_first = nCvor;
    } else {
        prev->next = nCvor;
    }

    return &nCvor->koeficijent;
}



std::string solve_exp(std::string base, std::string exp, bool baseInsidePara) {
    if ( isVariable(exp) ) {
        throw std::invalid_argument("Variable exponent not supported!");
    }

    std::vector<std::string> operands;
    Operators operators;
    splitExpression(base, operands, operators);

    if ( isVariable(base) && operands.size() == 1 ) {
        //ako je samo x na nesta, onda to ne diram nikako.
        //osim ako ima nesta uz x, onda to stepenujem prije nego sto vratim

        auto op = operands[0];
        //skini x

        auto var = op[op.size() - 1];
        removeVar(op);

        if ( baseInsidePara && op.size() > 0 ) {
            //imamo (Nx)^exp, pretvori u Kx^exp (gdje je K = N^exp)
            return solve_exp(op, exp, false) + var + "^" + exp;
        } else {
            //vrati samo baza^exp
            return op + var + "^" + exp;
        }
    }

    Polinom pa = Polinom::ParsirajPolinom(base);
    Polinom copy = pa;

    // stepenovanje rjesavamo sa visestrukim mnozenjem
    // npr a^5 pretvaramo u a * a * a * a * a

    int n = toNumber(exp);

    for ( int i = 1; i < n; ++ i ) {
        pa *= copy;
    }

    return pa.convertToString();
}

std::string solve_minus(std::string a, std::string b, bool aInsidePara) {
    auto pa = Polinom::ParsirajPolinom(a);
    auto pb = Polinom::ParsirajPolinom(b);

    auto pres = pa - pb;
    auto res = pres.convertToString();

    return res;
}

std::string solve_plus(std::string a, std::string b, bool aInsidePara) {
    auto pa = Polinom::ParsirajPolinom(a);
    auto pb = Polinom::ParsirajPolinom(b);

    auto pres = pa + pb;
    auto res = pres.convertToString();

    return res;
}

std::string solve_mul(std::string a, std::string b, bool aInsidePara) {
    if ( a == "0" || b == "0" ) {
        return "0";
    }

    auto pa = Polinom::ParsirajPolinom(a);
    auto pb = Polinom::ParsirajPolinom(b);

    auto pres = pa * pb;
    auto res = pres.convertToString();


    return res;
}


void splitExpression(std::string exp, std::vector<std::string>& operands, Operators& operators) {
    exp += " ";
    // std::cout<<"splitting -> "<<exp<<std::endl;

    std::string currExp;
    std::string paraContent;

    int paraCnt = 0;
    bool isMinus = false;

    if ( exp[0] == '-' ) {
        isMinus = true;
        exp.erase(exp.begin());
    }

    for ( int i = 0; i < exp.size(); ++ i ) {
        char c = exp[i];

        if ( isPara(c) ) {
            if ( c == '(' ) {
                if ( paraCnt > 0 ) {
                    paraContent += c;
                }
                ++paraCnt;
            } else {
                if ( paraCnt > 1 ) {
                    paraContent += c;
                }
                --paraCnt;
            }

            if ( paraCnt == 0 ) {
                Polinom paraPolinom = Polinom::ParsirajPolinom(paraContent);

                //ako je bio minus ispred zagrade, pomnozi sadrzaj zagrade sa -1 i prebaci
                // - u +

                if ( isMinus ) {
                    // ali samo ako zagrada nije uz operator ^, u tom slucaju mora prvo stepenovati pa onda negirati
                    if ( i < exp.size() - 1 && exp[i+1] == '^' ) {
                        //jeste zagrada uz ^, sredi to prvo...
                        i += 2; //preskoci zagradu i operator ^
                        std::string exponent;
                        while (i < exp.size() && isNum(exp[i]) ) {
                            exponent += exp[i];
                            ++i;
                        }
                        paraPolinom = Polinom::ParsirajPolinom(solve_exp(paraPolinom.convertToString(), exponent, true));
                    }
                    // pa onda nastavi normalno, pomnozi sa -1 zagradu i prebaci - u +
                    paraPolinom *= Polinom::ParsirajPolinom("-1");
                    if (operators.size()) {
                        operators[operators.size() - 1].op = '+';
                    }
                }
                operands.push_back(paraPolinom.convertToString());
                paraContent = "";
            }
        } else if (paraCnt != 0) {
            paraContent += c;
        } else if ( !isSeparator(c) ) {
            currExp += c;
        } else if ( isOperator(c) ) {
            bool nextOpIsMinus = c == '-';
            if (c == '-') {
                c = '+';
            }

            trim(currExp);
            if ( currExp.length() > 0 ) {
                // 5-3 se pretvara u 5 + (-3), ali
                // ne treba pretvoriti 5-3^2 u 5+(-3)^2
                if ( isMinus && (c != '^') ) {
                    currExp = '-' + currExp;
                } else if ( isMinus && (c == '^') ) {
                    //vrati zadnji operator u '-', da bude ipak 5-(3)^2
                    operators[operators.size() - 1].op = '-';
                }

                operands.push_back(currExp);
                currExp = "";
            }
            isMinus = nextOpIsMinus;

            // ne ubacuj kao operand ako je unarni.
            if ( operands.size() > 0 ) {
                operators.push_back(Operator(c, i));
            }
        } else {
            trim(currExp);
            if ( currExp.length() > 0 ) {
                if ( isMinus ) {
                    currExp = '-' + currExp;
                }

                operands.push_back(currExp);
                currExp = "";
            }
        }
    }

}
