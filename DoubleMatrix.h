/*
 *  DoubleMatrix.h
 *  
 *
 *  Created by Bin Wu on Oct/24/13.
 *
 */

//
// Matrix
//

class DoubleMatrix{
private:
	int _row,_col,*_indx;
	double **_v,**_lu,_d;
public:
	DoubleMatrix();
	DoubleMatrix(int m, int n);
	DoubleMatrix(const DoubleVector&);
	DoubleMatrix(const DoubleMatrix&);
	inline int row() const;
	inline int col() const;
	~DoubleMatrix();
	
	DoubleMatrix sub(int rowb,int rowe,int colb,int cole);
	DoubleVector rowVec(int rowb) const; 	//take out a row as a vector
	DoubleVector colVec(int colb) const;	//take out a col as a vector
	
	//Operators overloading
	inline double* operator[](const int i);
	inline const double* operator[](const int i) const;
	DoubleMatrix& operator=(const DoubleMatrix&);
	DoubleMatrix& operator=(const DoubleVector&);
	DoubleMatrix& operator=(const double&);
	DoubleMatrix operator+(const DoubleMatrix&) const;
	DoubleMatrix operator-(const DoubleMatrix&) const;
	DoubleMatrix operator*(const DoubleMatrix&) const;
	DoubleVector operator*(const DoubleVector&) const;	
	
	//Pointwise operation
	DoubleMatrix operator%(const DoubleMatrix&) const;
	DoubleMatrix operator/(const DoubleMatrix&) const;
	DoubleMatrix operator^(const double&) const;
	
	//LU Decomposition
	void LUdcmp();
	double det();
	DoubleVector operator|(const DoubleVector&) const;
	double integrate(const DoubleVector&) const;
	double integrate(const DoubleMatrix&) const;
	DoubleMatrix operator~() const;		//Take the inverse of the matrix
	DoubleMatrix operator!() const;		//Transport the matrix
};


/*****///Section 1: Definition

DoubleVector DoubleMatrix::rowVec(int rowb) const{
	DoubleVector rowV(_col);
	for(int i = 0;i<_col;i++) rowV[i] = _v[rowb][i];
	return rowV;
}

DoubleVector DoubleMatrix::colVec(int colb) const{
	DoubleVector colVec(_row);
	for(int i = 0;i<_row;i++) colVec[i] = _v[i][colb];
	return colVec;
}

DoubleMatrix::DoubleMatrix():_row(0),_col(0),_d(1.0),_v(NULL),_lu(NULL),_indx(NULL){}

DoubleMatrix::DoubleMatrix(int m, int n):_row(m),_col(n),_d(1.0),_lu(NULL),_indx(NULL),_v(m>0?new double*[m]:NULL){
	for(int i=0;i<_row;i++)
		_v[i]=_col>0?new double[_col]:NULL;
}

DoubleMatrix::DoubleMatrix(const DoubleMatrix& dm):_d(1.0),_indx(NULL){
	_row=dm.row();_col=dm.col();_lu=NULL;
	
	_v=_row>0?new double*[_row]:NULL;
	for(int i=0;i<_row;i++)
		_v[i]=_col>0?new double[_col]:NULL;
	
	for(int i=0;i<_row;i++)
		for(int j=0;j<_col;j++)
			_v[i][j]=dm[i][j];
}


DoubleMatrix::DoubleMatrix(const DoubleVector& v):_d(1.0),_indx(NULL){
	// Allocate new memory
	_row=v.size();_col=1;_lu=NULL;
		
	_v=_row>0?new double*[_row]:NULL;
	for(int i=0;i<_row;i++)_v[i]=new double[1];
	
	for(int i=0;i<_row;i++)
		_v[i][0]=v[i];
}

inline int DoubleMatrix::row() const{
	return _row;
}

inline int DoubleMatrix::col() const{
	return _col;
}

DoubleMatrix::~DoubleMatrix(){
	
	if(_lu!=NULL){
		for(int i=0;i<_row;i++) if(_lu[i]!=NULL) delete [] _lu[i];
		delete [] _lu;
	}
	
	if(_indx!=NULL) delete [] _indx;
	
	if(_v!=NULL){
		for(int i=0;i<_row;i++) if(_v[i]!=NULL) delete [] _v[i];
		delete [] _v;
	}
}

inline double* DoubleMatrix::operator[](const int i){
	if(i<0||i>=_row){
		cout << "***DoubleMatrix row out of bounds***\n";
		throw("DoubleMatrix row out of bounds");
	}
	return _v[i];
}

inline const double* DoubleMatrix::operator[](const int i) const{
	if(i<0||i>=_row){
		cout << "***DoubleMatrix row out of bounds***\n";
		throw("DoubleMatrix row out of bounds");
	}
	return _v[i];
}

DoubleMatrix DoubleMatrix::sub(int rowb,int rowe,int colb,int cole){
	if(rowe<_row&&cole<_col&&rowb>=0&&colb>=0){
		DoubleMatrix submatrix(rowe-rowb+1,cole-colb+1);
		int k=0;
		for(int i=rowb;i<=rowe;i++){
			int l=0;
			for(int j=colb;j<=cole;j++){
				submatrix[k][l]=_v[i][j];l++;
			}
			k++;
		}
		return submatrix;
	}
	else {
		cout << "***Out of bounds in subMatrix***\n";	
		throw("***Out of bounds in subMatrix***\n");
		return DoubleMatrix(0,0);
	}

}


/*****///Section 2: Operator overloading


DoubleMatrix& DoubleMatrix::operator=(const DoubleMatrix& v){
	if(this!=&v){
		if(_row!=v._row){
			// Free old memory
			if(_v!=NULL){
				for(int i=0;i<_row;i++) if(_v[i]!=NULL) delete [] _v[i];
				delete [] _v;
			}
			
			// Allocate new memory
			_row=v.row();_col=v.col();
			
			_v=_row>0?new double*[_row]:NULL;
			for(int i=0;i<_row;i++)
				_v[i]=_col>0?new double[_col]:NULL;
		}
		else if(_col!=v._col) {
			for(int i=0;i<_row;i++) if(_v[i]!=NULL) delete [] _v[i];
			_col=v._col;for(int i=0;i<_row;i++) _v[i]=_col>0?new double[_col]:NULL;
		}
		
		//cout << "***Start to release LU***\n";
		if(_lu!=NULL){
			for(int i=0;i<_row;i++) if(_lu[i]!=NULL) delete [] _lu[i];
			delete [] _lu;
		}
		
		if(_indx!=NULL) delete [] _indx;	
		
		_lu=NULL;_indx=NULL;

		//cout << "***End to release LU***\n";
		
		for(int i=0;i<_row;i++)
			for(int j=0;j<_col;j++)
				_v[i][j]=v[i][j];
	}
	return *this;
}

DoubleMatrix& DoubleMatrix::operator=(const DoubleVector& v){
	if(_row!=v.size()){
		// Free old memory
		if(_v!=NULL){
			for(int i=0;i<_row;i++) if(_v[i]!=NULL) delete [] _v[i];
			delete [] _v;
		}

		// Allocate new memory
		_row=v.size();_col=1;
			
		_v=_row>0?new double*[_row]:NULL;
		for(int i=0;i<_row;i++)
			_v[i]=new double[1];
	}
	else if(_col!=1) {
		for(int i=0;i<_row;i++) if(_v[i]!=NULL) delete [] _v[i];
		_col=1;for(int i=0;i<_row;i++) _v[i]=new double[_col];
	}		
	
	if(_lu!=NULL){
		for(int i=0;i<_row;i++) if(_lu[i]!=NULL) delete [] _lu[i];
		delete [] _lu;
	}
	
	if(_indx!=NULL) delete [] _indx;	

	_lu=NULL;_indx=NULL;
	
	for(int i=0;i<_row;i++)
		_v[i][0]=v[i];
	return *this;
}

DoubleMatrix& DoubleMatrix::operator=(const double& v){
	if(_v!=NULL){
		for (int i=0; i<_row; i++) for(int j=0;j<_col;j++) _v[i][j]=v;
	}else{
		cout << "***Give value to Empty Matrix***\n";	
		throw("***Give value to Empty Matrix***\n");		
	}
}


DoubleMatrix DoubleMatrix::operator*(const DoubleMatrix& v) const{
	DoubleMatrix nrv(_row,v.col());
	if (_col!= v.row()) {
		cout << "***Mutiplication of two Matrices with different dimensions***\n";
		throw("Mutiplication of two Matrices with different dimensions!\n");
	}
	else{
		for (int i=0; i<_row; i++){ 
			for (int j=0; j<v.col(); j++){
				nrv[i][j]=0;
				for(int k=0;k<_col;k++)
					nrv[i][j]+=(_v[i][k]*v[k][j]);
			}
		}
	}
	return nrv;
}

DoubleVector DoubleMatrix::operator*(const DoubleVector& v) const{
	DoubleVector nrv(_row);
	if (_col!= v.size()) {
		cout << "***Mutiplication of a Matrix and a Vector with different dimensions***\n";
		throw("Mutiplication of a Matrix and a Vector with different dimensions!\n");
	}
	else{
		for (int i=0; i<_row; i++){ 
				nrv[i]=0;
				for(int k=0;k<_col;k++)
					nrv[i]+=(_v[i][k]*v[k]);
		}
	}
	return nrv;
}

//***Note % mean pointwise multiplication

DoubleMatrix DoubleMatrix::operator%(const DoubleMatrix& v) const{
	DoubleMatrix nrv(_row,_col);
	if (_row!= v.row()||_col!= v.col()) {
		cout << "***PointWise multiplicate two Matrices with different dimensions***\n";
		throw("PointWise multiplicate two Matrices with different dimensions!\n");
	}
	else{
		for (int i=0; i<_row; i++) for (int j=0; j<_col; j++) nrv[i][j]=_v[i][j]*v[i][j];
	}
	return nrv;
}

//***Note % mean pointwise multiplication
DoubleMatrix DoubleMatrix::operator/(const DoubleMatrix& v) const{
	DoubleMatrix nrv(_row,_col);
	if (_row!= v.row()||_col!= v.col()) {
		cout << "***PointWise divide two Matrices with different dimensions***\n";
		throw("PointWise divide two Matrices with different dimensions!\n");
	}
	else{
		for (int i=0; i<_row; i++) for (int j=0; j<_col; j++) nrv[i][j]=_v[i][j]/v[i][j];
	}
	return nrv;
}

DoubleMatrix DoubleMatrix::operator+(const DoubleMatrix& v) const{
	DoubleMatrix nrv(_row,_col);
	if (_row!= v.row()||_col!= v.col()) {
		cout << "***Adding two Matrices with different dimensions***\n";
		throw("Adding two Matrices with different dimensions!\n");
	}
	else{
		for (int i=0; i<_row; i++) for (int j=0; j<_col; j++) nrv[i][j]=_v[i][j]+v[i][j];
	}
	return nrv;
}

DoubleMatrix DoubleMatrix::operator-(const DoubleMatrix& v) const{
	DoubleMatrix nrv(_row,_col);
	if (_row!= v.row()||_col!= v.col()) {
		cout << "***Substraction of two Matrices with different dimensions***\n";
		throw("Substraction of two Matrices with different dimensions!\n");
	}
	else{
		for (int i=0; i<_row; i++) for (int j=0; j<_col; j++) nrv[i][j]=_v[i][j]-v[i][j];
	}
	return nrv;
}

/*****///Section 3: Linear Algebra

// LU Decomposition from numerical recipies 3rd: M = LU
void DoubleMatrix::LUdcmp(){
	if(_row==_col){
		if(_lu==NULL){
			
			const double TINY=1.0e-40;
			int i,imax,j,k,n=_row;
			double big,temp;
			DoubleVector vv(n);
			
			_lu=_row>0?new double*[_row]:NULL;
			for(i=0;i<_row;i++)
				_lu[i]=_col>0?new double[_col]:NULL;
			
			for(i=0;i<_row;i++)
				for(j=0;j<_col;j++)
					_lu[i][j]=_v[i][j];	
			
			_indx=new int[_row];

			_d=1.0;
			for (i=0;i<n;i++) {
				big=0.0;
				for (j=0;j<n;j++)
					if ((temp=abs(_lu[i][j])) > big) big=temp;
				if (big == 0.0){
					cout << "***Singular matrix in LUdcmp***\n";
					throw("Singular matrix in LUdcmp");
				}
				vv[i]=1.0/big;
			}
			for (k=0;k<n;k++) {
				big=0.0;
				for (i=k;i<n;i++) {
					temp=vv[i]*abs(_lu[i][k]);
					if (temp > big) {
						big=temp;
						imax=i;
					}
				}
				if (k != imax) {
					for (j=0;j<n;j++) {
						temp=_lu[imax][j];
						_lu[imax][j]=_lu[k][j];
						_lu[k][j]=temp;
					}
					_d = -_d;
					vv[imax]=vv[k];
				}
				_indx[k]=imax;
				if (_lu[k][k] == 0.0) _lu[k][k]=TINY;
				for (i=k+1;i<n;i++) {
					temp=_lu[i][k] /= _lu[k][k];
					for (j=k+1;j<n;j++)
						_lu[i][j] -= temp*_lu[k][j];
				}
			}
		}
	}
}

//Solve the linear equation x = M^(-1)*b = U^-1*L^-1*b 
DoubleVector DoubleMatrix::operator|(const DoubleVector& b) const
{
	int i,ii=0,ip,j;
	double sum;int n=_row;
	DoubleVector x(_row);
	
	if (b.size() != _row){
		cout << "***LUdcmp::solve bad sizes***\n";
		throw("LUdcmp::solve bad sizes");
	}
	
	for (i=0;i<n;i++) x[i] = b[i];
	
	//L^{-1}*b
	for (i=0;i<n;i++) {
		ip=_indx[i];
		sum=x[ip];
		x[ip]=x[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= _lu[i][j]*x[j];
		else if (sum != 0.0)
			ii=i+1;
		x[i]=sum;
	}

	//U^{-1}*(L^{-1}*b)
	for (i=n-1;i>=0;i--) {
		sum=x[i];
		for (j=i+1;j<n;j++) sum -= _lu[i][j]*x[j];
		x[i]=sum/_lu[i][i];
	}
	
	return x;
}

//Integration of b
double DoubleMatrix::integrate(const DoubleVector& b) const
{
	int i,ii=0,ip,j;
	double sum;int n=_row;
	DoubleVector x=b;
	
	if (b.size() != _row){
		cout << "***LUdcmp::solve bad sizes***\n";
		throw("LUdcmp::solve bad sizes");
	}
	
	//L^{-1}*b
	for (i=0;i<n;i++) {
		ip=_indx[i];
		sum=x[ip];
		x[ip]=x[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= _lu[i][j]*x[j];
		else if (sum != 0.0)
			ii=i+1;
		x[i]=sum;
	}
	
	i=n-1;
	sum = x[i]/_lu[i][i];
	
	return sum;
}

//Integration of Matrix b: here I will have to understand why const for rowVec works???
double DoubleMatrix::integrate(const DoubleMatrix& b) const
{
	int i; DoubleVector x(b.row());

	for(i=0;i<b.row();i++){
	 	x[i]=integrate(b.rowVec(i));
	}
	
	return integrate(x);
}

double DoubleMatrix::det()
{
	double dd = _d;
	for (int i=0;i<_row;i++) dd *= _lu[i][i];
	return dd;
}

DoubleMatrix DoubleMatrix::operator!() const{	//Transport the matrix
	DoubleMatrix m(_row,_col);int i,j;
	for(i=0;i<_row;i++)
		for(j=0;j<_col;j++)
			m[i][j] = _v[j][i];

	return m;
}


//Inverse of M
DoubleMatrix DoubleMatrix::operator~() const{
	DoubleVector x(_row),y;
	DoubleMatrix invm(_row,_col);
	
	for(int i=0;i<_col;i++){
		for(int j=0;j<_row;j++) x[j]=0;x[i]=1;
		y=(this->operator|(x));
		for(int j=0;j<_row;j++) invm[j][i]=y[j];
	}
	
	return invm;
}



