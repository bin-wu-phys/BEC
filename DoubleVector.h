/*
 *  DoubleVector.h
 *  
 *
 *  Created by Bin Wu on 5/25/12.
 *  Copyright 2012 Universitaet Bielefeld. All rights reserved.
 *
 */
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
//#define _MATLAB

//
//
// Vector
//
//
#define PI 3.141592653589793//23846 26433 83279 50288 41971 69399

class DoubleVector{
private:
	int _n;
	double *_v;

public:
	//Class definition
	DoubleVector();
	DoubleVector(int);
	DoubleVector(double v);
	DoubleVector(int n,double *v);
	DoubleVector(const DoubleVector& v);

	~DoubleVector();
	int size() const;

	DoubleVector sub(int rowb,int rowe); //Get subVector
	
	//Operators overloading
	double& operator[](int);
	const double& operator[](int) const;
	DoubleVector& operator=(const DoubleVector&);
	DoubleVector& operator=(const double);
	DoubleVector operator+(const DoubleVector&) const;
	DoubleVector operator-(const DoubleVector&) const;
	DoubleVector operator*(const DoubleVector&) const;
	DoubleVector operator/(const DoubleVector&) const;
	DoubleVector operator^(const double&) const;
	
	//Data analysis
	double sum();			//get  sum of the elements of the vector
	double last();			//get  the last element of the vector
	DoubleVector every(int sub);	//get  the sub vector take _n[i mod sub]
};

DoubleVector::DoubleVector():_n(0),_v(NULL){}

DoubleVector::DoubleVector(int n){
	_v=new double[n];_n=n;
}

DoubleVector::DoubleVector(double v){
	_v=new double[1];_n=1;_v[0]=v;
}

DoubleVector::DoubleVector(int n,double *v){
	_v=new double[n];_n=n;for(int i=0;i<n;i++) _v[i]=v[i];
}

DoubleVector::DoubleVector(const DoubleVector& v){
	_n=v._n;
	_v=_n>0?new double[_n]:NULL;
	for(int i=0;i<_n;i++) _v[i]=v[i];
}

int DoubleVector::size() const{//const is also needed!
	return _n;
}

DoubleVector::~DoubleVector(){
	if(_v!=NULL) delete[] _v;
}

double& DoubleVector::operator[](int i){
	if(i<0||i>_n){
		cout << "Out of Bounds of DoubleVector\n";
		throw string("Out of Bounds of DoubleVector\n");
	}
	return _v[i];
}

const double& DoubleVector::operator[](int i) const{
	if(i<0||i>_n){
		cout << "Out of Bounds of DoubleVector\n";
		throw string("Out of Bounds of DoubleVector\n");
	}
	return _v[i];
}

DoubleVector& DoubleVector::operator=(const DoubleVector& v){
	if(this!=&v){
		if(_n!=v._n){
			if(_v!=NULL) delete[] _v;
			_n=v._n;
			_v=_n>0?new double[_n]:NULL;
		}
		for(int i=0;i<_n;i++)
			_v[i]=v[i];
	}
	return *this;
}

DoubleVector& DoubleVector::operator=(const double v){
	if(_n!=1){
		if(_v!=NULL) delete[] _v;
		_n=1;
		_v=_n>0?new double[_n]:NULL;
	}
		_v[0]=v;
	return *this;
}

DoubleVector DoubleVector::operator+(const DoubleVector& v) const{
	DoubleVector nrv(_n);
	if (_n!= v.size()) {
		cout << "***Adding two vectors with different dimensions!\n***";
		throw("Adding two vectors with different dimensions!\n");
	}
	else{
		for (int i=0; i<_n; i++) nrv[i]=_v[i]+v[i];
	}
	return nrv;
}

DoubleVector DoubleVector::operator-(const DoubleVector& v) const{
	DoubleVector nrv(_n);
	if (_n!= v.size()) {
		cout << "***Substraction of two vectors with different dimensions!\n***";
		throw("Substraction of two vectors with different dimensions!\n");
	}
	else{
		for (int i=0; i<_n; i++) nrv[i]=_v[i]-v[i];
	}
	return nrv;
}

DoubleVector DoubleVector::operator*(const DoubleVector& v) const{
	DoubleVector nrv(_n);
	if (_n!= v.size()) {
		cout << "***Multiplication of two vectors with different dimensions!\n***";
		throw("Multiplication of two vectors with different dimensions!\n");
	}
	else{
		for (int i=0; i<_n; i++) nrv[i]=_v[i]*v[i];
	}
	return nrv;
}

DoubleVector DoubleVector::operator/(const DoubleVector& v) const{
	DoubleVector nrv(_n);
	if (_n!= v.size()) {
		cout << "***Division of two vectors with different dimensions!\n***";
		throw("Division of two vectors with different dimensions!\n");
	}
	else{
		for (int i=0; i<_n; i++) nrv[i]=_v[i]/v[i];
	}
	return nrv;
}

DoubleVector DoubleVector::operator^(const double& v) const{
	DoubleVector nrv(_n);
	for (int i=0; i<_n; i++) nrv[i]=pow(_v[i],v);
	return nrv;
}


DoubleVector DoubleVector::sub(int rowb,int rowe){
	DoubleVector subVector(rowe-rowb+1);
	if(rowe<_n&&rowb>=0){
		int k=0;
		for(int i=rowb;i<=rowe;i++){
			subVector[k]=_v[i];k++;
		}
	}else{
		cout << "***Out of bounds in subVector***\n";
		throw("***Out of bounds in subVector***\n");
	}
	return subVector;	
}

//Data analysis
double DoubleVector::sum(){
	double sum=0;
	for(int i=0;i<_n;i++) sum+=_v[i];
	return sum;
}

double DoubleVector::last(){
	return _v[_n-1];
}
