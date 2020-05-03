#pragma once
#include <cmath>
#include <complex>

namespace solver{

    
    
    class RealVariable{
        public:
            double a;
            double b;
            double c;


        //ax^2+bx+c
            RealVariable();
        RealVariable(double a, double b, double c);
        RealVariable(const RealVariable& copy);

        //+
        friend RealVariable& operator+(double y,RealVariable& x);
        friend RealVariable& operator+(RealVariable& x, double y);
        friend RealVariable& operator+(RealVariable& x,RealVariable& y);
        //*
        friend RealVariable& operator*(double y,RealVariable& x);
        friend RealVariable& operator*(RealVariable& x, double y);
        friend RealVariable& operator*(RealVariable& x,RealVariable& y);
        ///
        friend RealVariable& operator/(RealVariable& x,double y);
        friend RealVariable& operator/(double y,RealVariable& x);
        friend RealVariable& operator/(RealVariable& x,RealVariable& y);
        //-
        friend RealVariable& operator-(RealVariable& x, double y);
        friend RealVariable& operator-(RealVariable& x,RealVariable& y);
        friend RealVariable& operator-(double y,RealVariable& x);
        friend RealVariable& operator-(RealVariable& x);
        //^
        friend RealVariable& operator^(RealVariable& x, double y);
        //==
        friend RealVariable& operator==(RealVariable& x, double y);
        friend RealVariable& operator==(RealVariable& x,RealVariable& y);
        friend RealVariable& operator==(double y,RealVariable& x);
        
    
    };
    double solve(RealVariable& x);
    
    class ComplexVariable{
        public:
        std::complex<double> a,b,c;
        ComplexVariable()=default;

        friend ComplexVariable& operator+(int y,ComplexVariable& x);
        friend ComplexVariable& operator+(ComplexVariable& x,int y);
        friend ComplexVariable& operator+(ComplexVariable& x,ComplexVariable& y);
        friend ComplexVariable& operator+(std::complex<double> y,ComplexVariable& x);
        friend ComplexVariable& operator+(ComplexVariable& x,std::complex<double> y);

        //*
        friend ComplexVariable& operator*(int y,ComplexVariable& x);
        friend ComplexVariable& operator*(ComplexVariable& x, int y);
        friend ComplexVariable& operator*(ComplexVariable& x,ComplexVariable& y);
        ///
        friend ComplexVariable& operator/(ComplexVariable& x, double y);
        friend ComplexVariable& operator/(double y,ComplexVariable& x);
        friend ComplexVariable& operator/(ComplexVariable& x,ComplexVariable& y);
        //-
        friend ComplexVariable& operator-(ComplexVariable& x,int y);
        friend ComplexVariable& operator-(ComplexVariable& x,ComplexVariable& y);
        friend ComplexVariable& operator-(int y,ComplexVariable& x);
        //^
        friend ComplexVariable& operator^(ComplexVariable& x,int y);
        //==
        friend ComplexVariable& operator==(ComplexVariable& x, int y); 
        friend ComplexVariable& operator==(ComplexVariable& x,ComplexVariable& y);
        friend ComplexVariable& operator==(int y,ComplexVariable& x);
        friend ComplexVariable& operator==(ComplexVariable& x,std::complex<double> y);
        friend ComplexVariable& operator==(std::complex<double> y,ComplexVariable& x);

    
    };
    std::complex<double> solve(ComplexVariable& x);

}