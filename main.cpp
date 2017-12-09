#include <iostream>
#include "Polinom.h"


void testiraj(const std::string& expr, const std::string& res) {
    Polinom p = Polinom::ParsirajPolinom(expr);

    if ( p.convertToString() != res ) {
        std::cout<<"Dobila '"<<p.convertToString()<<"', a trebala '"<<res<<"'"<<std::endl;
        throw std::invalid_argument("Neispravan rezultat!");
    } else {
        std::cout<<expr<<" = "<<res<<" // OK "<<std::endl;
    }
}

int main() {
    /*//testovi bez nepoznatih
    testiraj("2 + 5 - 2^3 - 2", "-3");
    testiraj("2 + (5 - 2)^3 - 2", "27");
    testiraj("2 + (5 - (2-4+7))^3 - 2^6", "-62");


    //testovi sa nepoznatom
    testiraj("(0+3x)^2 + 3x^2", "12x^2");
    testiraj("(2+3x)^3 + 3x^2", "8 + 36x + 57x^2 + 27x^3");

    testiraj("2x + (5 - 3)^3 - 2", "6 + 2x");
    testiraj("2x + (5 - 3)^3 - 2 - 5x^4", "6 + 2x - 5x^4");
    testiraj("-x-1", "-1 - 1x");

    testiraj("(x+(x-1)^2)^3", "1 - 3x + 6x^2 - 7x^3 + 6x^4 - 3x^5 + 1x^6");
    testiraj("-(x+(x-1)^2)^3", "-1 + 3x - 6x^2 + 7x^3 - 6x^4 + 3x^5 - 1x^6");

    testiraj("((2+(x+9)^2)*(x+1))", "83 + 101x + 19x^2 + 1x^3");


    testiraj("(3x-(-x-1)^2)*(-2x+3)", "-3 + 5x - 5x^2 + 2x^3");

    Polinom p1 = Polinom::ParsirajPolinom("(x+(x-1)^2)^3");
    Polinom p2 = Polinom::ParsirajPolinom("(x+(x-1)^2)^2");

    std::cout<<std::endl;
    std::cout<<"("<<p1<<") * ("<<p2<<") = "<<(p1 * p2)<<std::endl;
    std::cout<<"("<<p1<<") + ("<<p2<<") = "<<(p1 + p2)<<std::endl;
    std::cout<<"("<<p1<<") - ("<<p2<<") = "<<(p1 - p2)<<std::endl;

*/
    // proba unosa
    try {
        std::cout<<"Unesite izraz:";
        Polinom p3;
        std::cin >> p3;
        std::cout << p3 << std::endl;
    } catch ( const std::exception& ex ) {
        std::cout<<ex.what();
    }



    return 0;
}

