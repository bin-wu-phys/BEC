/*
 *  This file includes all the top-level operator overloading of the DoubleVector class
 *  
 *
 *  Created by Bin Wu on May/25/12.
 *  Modified on Jan 3, 2014.
 */
//
//Top-level functions
//
#define INF     1e20

//Output
ostream& operator<<(ostream& out,const DoubleVector& v){
	out.precision(std::numeric_limits<double>::digits10);
	out << setiosflags(ios::scientific | ios::showpoint);
	if(v.size()>0){
#ifdef _MATLAB
		int nm1=v.size()-1;
		out <<" [ ";
		for(int i=0;i<nm1;i++) out << v[i] <<"; ";
		out << v[nm1] << " ] ";
#else
	#ifdef _MATHEMATICA
		int last = v.size()-1;
		for(int i=0;i<v.size();i++){
			if(i!=last) out << v[i] << "   ";
			else  out << v[i];
		}

	#else
		out << "\n";
		for(int i=0;i<v.size();i++){
			out << v[i] << "\n";
		}
		out << "\n";
	#endif
#endif
	}else {
		out << "NULL";
	}
	
	return out;
}

DoubleVector operator&(const DoubleVector& v1,const DoubleVector& v2){
	
	if(v1.size()>0&&v2.size()>0){
		DoubleVector nrv(v1.size()+v2.size());
		for(int j=0;j<v1.size();j++) nrv[j]=v1[j];
		for(int k=0;k<v2.size();k++) nrv[v1.size()+k]=v2[k];
		return nrv;
	}
	else {
		if(v1.size()<=0) return v2;
		else if(v2.size()<=0) return v1;
		else{
			return DoubleVector(0);
		}
	}	
}


DoubleVector operator&(const double& v1,const DoubleVector& v2){
	
	if(v2.size()>0){
		DoubleVector nrv(1+v2.size());
		nrv[0]=v1;
		for(int k=0;k<v2.size();k++) nrv[1+k]=v2[k];
		return nrv;
	}
	else {
		return DoubleVector(v1);
	}	
}

DoubleVector operator&(const DoubleVector& v1,const double& v2){
	
	if(v1.size()>0){
		DoubleVector nrv(1+v1.size());
		for(int k=0;k<v1.size();k++) nrv[k]=v1[k];
		nrv[v1.size()]=v2;
		return nrv;
	}
	else {
		return DoubleVector(v2);
	}	
}

// Top-level operator overloading

DoubleVector operator+(const double& dv,const DoubleVector& v){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=dv+v[i];
	return nrv;	
}

DoubleVector operator+(const DoubleVector& v, const double& dv){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=dv+v[i];
	return nrv;	
}

DoubleVector operator-(const double& dv,const DoubleVector& v){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=dv-v[i];
	return nrv;	
}

DoubleVector operator-(const DoubleVector& v){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=-v[i];
	return nrv;	
}

DoubleVector operator-(const DoubleVector& v, const double& dv){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=v[i]-dv;
	return nrv;	
}

DoubleVector operator*(const double& dv,const DoubleVector& v){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=dv*v[i];
	return nrv;	
}

DoubleVector operator/(const double& dv,const DoubleVector& v){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=dv/v[i];
	return nrv;	
}

DoubleVector operator/(const DoubleVector& v, const double& dv){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=v[i]/dv;
	return nrv;	
}

//Scalar product
double operator%(const DoubleVector& dv,const DoubleVector& v){
	double sum=0;
	for (int i=0; i<v.size(); i++) sum+=dv[i]*v[i];
	return sum;	
}


// Mathematical functions overloading

DoubleVector sin(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){
		if(v[i]==PI)
			retV[i]=0;
		else
			retV[i]=sin(v[i]);
	}
	return retV;
}

DoubleVector exp(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=exp(v[i]);
	}
	return retV;
}

DoubleVector log(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){
		if(v[i]==0)
			retV[i]=-INF;
		else
			retV[i]=log(abs(v[i]));
	}
	return retV;
}


DoubleVector pow(const DoubleVector& v,const double &ex){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=pow(v[i],ex);
	}
	return retV;
}

DoubleVector cos(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		if(v[i]==0.5*PI)
			retV[i]=0;
		else
			retV[i]=cos(v[i]);
	}
	return retV;
}

DoubleVector tan(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=tan(v[i]);
	}
	return retV;
}

DoubleVector cot(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		if(v[i]==0.5*PI)
			retV[i]=0;
		else
			retV[i]=1/tan(v[i]);
	}
	return retV;
}

DoubleVector asin(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=asin(v[i]);
	}
	return retV;
}

DoubleVector acos(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=acos(v[i]);
	}
	return retV;
}



