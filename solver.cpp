#include "solver.hpp"
using namespace std;
using namespace solver;
#include <math.h>  



double solver::solve(RealVariable& x){

    if (x.a == 0 && x.b == 0)
        throw invalid_argument("");
    
    if (x.a != 0)
    {
        double d = (x.b * x.b) + (-4 * x.a * x.c);
        double solution = 0;

        if (d < 0)
            throw invalid_argument("");
        else
        {
            double sqrtD = sqrt(d);
            solution = (-x.b);
            solution += sqrtD;
            solution = solution / (2 * x.a);
        }

        return solution;
    }
    else
    {
        return (-x.c) / x.b;
    }
    
}

solver::RealVariable::RealVariable() {
    this->a = 0;
    this->b = 1;
    this->c = 0;
}
solver::RealVariable::RealVariable(double a, double b, double c)
{
    this->a = a;
    this->b = b;
    this->c = c;
}

RealVariable::RealVariable(const RealVariable& copy)
{
    this->a = copy.a;
    this->b = copy.b;
    this->c = copy.c;
}

RealVariable& solver::operator+(RealVariable& x,RealVariable& y){
    RealVariable* z = new RealVariable(x);
    z->a = x.a + y.a;
    z->b= x.b + y.b;
    z->c = x.c + y.c;
    return *z;
}

RealVariable& solver::operator+(RealVariable& x, double y){
    RealVariable* z = new RealVariable(x);
    z->c = x.c + y;
    return *z;

}
RealVariable& solver::operator+(double y,RealVariable& x){

    RealVariable* temp = new RealVariable(x);
    RealVariable z  =*temp + y;
    *temp = z;
    return *temp;

}

RealVariable& solver::operator*(double x,RealVariable& y){
    RealVariable* temp =new RealVariable(0, 0, x);
    RealVariable z = *temp * y;
    *temp = z;
    return *temp;
}

RealVariable& solver::operator*(RealVariable& x,RealVariable& y){

    RealVariable* z = new RealVariable(x);
    z->a = x.a * y.c + x.b * y.b + x.c * y.a;
    z->b = x.b * y.c + x.c * y.b;
    z->c = x.c * y.c;
    return *z;

}

RealVariable& solver::operator*(RealVariable& y, double x){
    RealVariable temp = RealVariable(0, 0, x);
    return y * temp;
}
RealVariable& solver::operator/(RealVariable& x,double y){
    if (y == 0)
        throw invalid_argument("");
    
    RealVariable temp = RealVariable(0, 0, 1 / y);
    return x * temp;

}
RealVariable& solver::operator/(RealVariable& x,RealVariable& y){
    
    return x;

}
RealVariable& solver::operator/(double x,RealVariable& y){
        
    return y;

}
RealVariable& solver::operator-(RealVariable& x, double y){
    RealVariable temp = RealVariable(0, 0, -y);
    return x + temp;

}

RealVariable& solver::operator-(RealVariable& x)
{
    return -1 * x;
}

RealVariable& solver::operator-(double x,RealVariable& y){
        
    RealVariable* temp = new RealVariable(y);
    *temp = -(*temp);
    RealVariable* temp2 = new RealVariable(0, 0, 0);
    *temp2 = x + y;
    delete temp;
    return *temp2;

}
RealVariable& solver::operator-(RealVariable& x,RealVariable& y){
       
    RealVariable temp = -y;
    return x + temp;
}
        
RealVariable& solver::operator^(RealVariable& x, double y){
    
    if (y > 2 || y < 0)
        throw invalid_argument("");
    RealVariable* temp = new RealVariable(x);
        if (y == 0)
        {
            temp->a = 0;
            temp->b = 0;
            temp->c = 1;
        }

        for (int i = 1; i < y; i++)
            *temp = x * x;

        return *temp;
} 
RealVariable& solver::operator==(RealVariable& x, RealVariable& y) {
        
    RealVariable* minusY = new RealVariable(y);
    *minusY = -y;
    *minusY = *minusY + x;
    return *minusY;
}
RealVariable& solver::operator==(double x,RealVariable& y) {
        
    RealVariable* minusX = new RealVariable(0,0,-x);
    RealVariable* result = &(y + *minusX);
    delete minusX;
    return *result;

}
RealVariable& solver::operator==(RealVariable& x, double y) {

    return y == x;
        
}

//******************

std::complex<double> solver::solve(ComplexVariable& x){
    return 0;

}

ComplexVariable& solver::operator+(ComplexVariable& x,ComplexVariable& y){
    return x;
}

ComplexVariable& solver::operator+(ComplexVariable& x,int y){
    
    return x;

}
ComplexVariable& solver::operator+(int y,ComplexVariable& x){

    return x;

}
ComplexVariable& solver::operator+(std::complex<double> y,ComplexVariable& x){
    return x;
}


ComplexVariable& solver::operator+(ComplexVariable& x,std::complex<double> y){
    return x;
}

ComplexVariable& solver::operator*(int x,ComplexVariable& y){
    return y;
}

ComplexVariable& solver::operator*(ComplexVariable& x,ComplexVariable& y){
    return y;
}

ComplexVariable& solver::operator*(ComplexVariable& y,int x){
    return y;
}
ComplexVariable& solver::operator/(ComplexVariable& x,double y){
    return x;

}
ComplexVariable& solver::operator/(ComplexVariable& x,ComplexVariable& y){
        
    return x;

}
ComplexVariable& solver::operator/(double x,ComplexVariable& y){
        
    return y;

    }
ComplexVariable& solver::operator-(ComplexVariable& x, int y){
    return x;

    }
ComplexVariable& solver::operator-(int x,ComplexVariable& y){
        
    return y;
}
ComplexVariable& solver::operator-(ComplexVariable& x,ComplexVariable& y){
        
    return x;
}
        
ComplexVariable& solver::operator^(ComplexVariable& x,int y){
    return x;

} 
ComplexVariable& solver::operator==(ComplexVariable& x, ComplexVariable& y) {
        
    return x;
}
ComplexVariable& solver::operator==(int x,ComplexVariable& y) {
        
    return y;
}
ComplexVariable& solver::operator==(ComplexVariable& x, int y) {

    return x;
        
}
ComplexVariable& solver::operator==(ComplexVariable& x,std::complex<double>){

    return x;
}

ComplexVariable& solver::operator==(std::complex<double> y,ComplexVariable& x){
    return x;
}


