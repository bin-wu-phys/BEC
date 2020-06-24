/*
 *  This file includes all the top-level operator overloading of the DoubleMatrix class
 *  
 *
 *  Created by Bin Wu on Oct/24/13.
 *
 */

// Top-level operator overloading

ostream& operator<<(ostream& out,const DoubleMatrix& m){
        out.precision(std::numeric_limits<double>::digits10);
        out << setiosflags(ios::scientific |
ios::showpoint);
	if(m.row()>0&&m.col()>0){
#ifdef _MATLAB
		int nrow=m.row()-1, ncol=m.col()-1;
		out <<" [";
		for(int i=0;i<nrow;i++){
			//cout << "{ ";
			for(int j=0;j<ncol;j++)
				out << m[i][j] <<", ";
			out << m[i][ncol] << " ; ";
		}
		//cout << "{ ";
		for(int j=0;j<ncol;j++)
			out << m[nrow][j] <<", ";
		out << m[nrow][ncol] << " ] ";
#else
	#ifdef _MATHEMATICA
		int collast = m.col()-1,rowlast = m.row()-1;
		for(int i=0;i<m.row();i++){
			for(int j=0;j<m.col();j++){
				if(j!=collast)
					out << m[i][j] << "   ";
				else
					out << m[i][j];	
			}
			if(i<rowlast) out << "\n";		
		}
	 
	#else
		out<<"\n";
		for(int i=0;i<m.row();i++){
			for(int j=0;j<m.col();j++){
				out << m[i][j] << "   ";
			}
			out << "\n";		
		}
		out << "\n";
	#endif
#endif
	}else {
		out << "NULL";
	}
	return out;
}

DoubleMatrix operator*(const double& dv,const DoubleMatrix& v){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for(int j=0;j<v.col();j++) nrv[i][j]=dv*v[i][j];
	return nrv;	
}

DoubleMatrix operator*(const DoubleMatrix& v, const double& dv){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for(int j=0;j<v.col();j++) nrv[i][j]=dv*v[i][j];
	return nrv;	
}

DoubleMatrix operator*(const DoubleVector& dv,const DoubleMatrix& v){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for(int j=0;j<v.col();j++) nrv[i][j]=dv[i]*v[i][j];
	return nrv;	
}

DoubleMatrix operator+(const DoubleMatrix& m,const DoubleVector& v){
	DoubleMatrix nrv=m;
	if(m.row()==m.col()&&m.row()==v.size()){
		for (int i=0; i<m.row(); i++) nrv[i][i]+=v[i];
	}else {
		cout << "***Addition of Matrix and Vector with different dimensions***\n";
		throw("***Addition of Matrix and Vector with different dimensions***\n");
	}

	return nrv;	
}

DoubleMatrix operator+(const double& dv,const DoubleMatrix& v){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for (int j=0; j<v.col(); j++) nrv[i][j]=dv+v[i][j];
	return nrv;
}

DoubleMatrix operator+(const DoubleMatrix& v, const double& dv){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for (int j=0; j<v.col(); j++) nrv[i][j]=dv+v[i][j];
	return nrv;
}

DoubleMatrix operator+(const DoubleVector& v,const DoubleMatrix& m){
	DoubleMatrix nrv=m;
	if(m.row()==m.col()&&m.row()==v.size()){
		for (int i=0; i<m.row(); i++) nrv[i][i]+=v[i];
	}else {
		cout << "***Addition of Matrix and Vector with different dimensions***\n";
		throw("***Addition of Matrix and Vector with different dimensions***\n");
	}
	
	return nrv;	
}


DoubleMatrix operator-(const double& dv, const DoubleMatrix& v){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for (int j=0; j<v.col(); j++) nrv[i][j]=dv-v[i][j];
	return nrv;
}

DoubleMatrix operator-(const DoubleMatrix& v, const double& dv){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for (int j=0; j<v.col(); j++) nrv[i][j]=v[i][j]-dv;
	return nrv;
}

DoubleMatrix operator-(const DoubleMatrix& m){
	DoubleMatrix nrv(m.row(),m.col());
	for (int i=0; i<m.row(); i++) for(int j=0;j<m.col();j++) nrv[i][j]=-m[i][j];
	return nrv;	
}

DoubleMatrix operator-(const DoubleMatrix& m,const DoubleVector& v){
	DoubleMatrix nrv=m;
	if(m.row()==m.col()&&m.row()==v.size()){
		for (int i=0; i<m.row(); i++) nrv[i][i]-=v[i];
	}else {
		cout << "***Subtraction of Matrix and Vector with different dimensions***\n";
		throw("***Subtraction of Matrix and Vector with different dimensions***\n");
	}
	
	return nrv;	
}

DoubleMatrix operator-(const DoubleVector& v,const DoubleMatrix& m){
	DoubleMatrix nrv=(-m);
	if(m.row()==m.col()&&m.row()==v.size()){
		for (int i=0; i<m.row(); i++) nrv[i][i]+=v[i];
	}else {
		cout << "***Subtraction of Matrix and Vector with different dimensions***\n";
		throw("***Subtraction of Matrix and Vector with different dimensions***\n");
	}
	
	return nrv;	
}

DoubleMatrix operator/(const DoubleMatrix& v,const double& dv){
	DoubleMatrix nrv(v.row(),v.col());
	if(dv!=0){
		for (int i=0; i<v.row(); i++) for(int j=0;j<v.col();j++) nrv[i][j]=v[i][j]/dv;
	}else {
		cout << "***ZERO divied by Matrix***\n";
		throw("***ZERO divied by Matrix***\n");
	}

	return nrv;	
}


DoubleMatrix operator/(const double& dv, const DoubleMatrix& v){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for(int j=0;j<v.col();j++){ 
		if(v[i][j]!=0){
			nrv[i][j]=dv/v[i][j];
			}else {
				cout << "***ZERO Element of Matrix divied by a double num***\n";
				throw("***ZERO Element of Matrix divied by a double num***\n");
			}
	}
	return nrv;	
}

DoubleMatrix operator&(const DoubleMatrix& v1,const DoubleMatrix& v2){
	
	if(v1.row()==v2.row()){
		DoubleMatrix nrv(v1.row(),v1.col()+v2.col());
		for(int i=0; i<v1.row(); i++){
			for(int j=0;j<v1.col();j++) 
				nrv[i][j]=v1[i][j];
			for(int k=0;k<v2.col();k++)
				nrv[i][v1.col()+k]=v2[i][k];
		}
		return nrv;
	}
	else {
		if(v1.row()<=0) return v2;
		else if(v2.row()<=0) return v1;
		else{
			cout << "***Trying to combine Matrices with different dimensions***\n";
			return DoubleMatrix(0,0);
		}
	}	
}

#ifdef _MATLAB
void OutputList(const DoubleMatrix& m){
	for(int i=0;i<m.row();i++){
		for(int j=0;j<m.col();j++){
			cout << m[i][j] << "   ";
		}
		cout << "\n";		
	}
}
#endif
