using namespace std;
const double Ee=2.71828182845904523536; //e
const double PI=3.14159265358979323846; //pi
double det(const int n,const double number[ ]); //计算实系数方阵的行列式
double tr(const int n,const double number[ ]); //计算矩阵的迹
double rank(const int r,const int c,const double number[ ]); //求矩阵秩
double fact(const double x,bool *matherror); //计算阶乘
double factdouble(const double x,bool *matherror); //计算双阶乘
double count(double x,const double y,const int z,bool *matherror,const bool degreemode); //对某种具体运算实行计算，调用了多个计算的函数
inline bool cmp(const double x,const double y,const int z); //比较大小
double gcd(double x,double y,bool *matherror); //计算浮点型最大公约数
int gcdn(int x,int y); //计算整型最大公约数
double lcm(const double x,const double y,bool *matherror); //计算最小公倍数
double ran(const double x); //找随机数
inline double minnum(const double x,const double y); //计算两数最小值
inline double maxnum(const double x,const double y); //计算两数最大值
double P(const double x,const double y,bool *matherror); //计算排列数
double C(const double x,const double y,bool *matherror); //计算组合数
double gauss(const double x); //高斯函数（向下取整）
inline bool isint(const double x); //判断一个数是否整数
double __pow(const double x,double y,bool *matherror); //幂函数（C++原配的pow函数的改进版）
inline double loginot(const double x); //逻辑非运算
inline double logiand(const double x,const double y); //逻辑与运算
inline double logior(const double x,const double y); //逻辑或运算
inline double logixor(const double x,const double y); //逻辑异或运算
double fib(double x,bool *matherror); //斐波那契数列通项
inline short sgn(const double x); //符号函数
double round(double x,const double y); //四舍五入
double factornumber(const double x); //求约数的总数
double planeangle(const double x,const double y); //计算极角
int daysbetween(int year1,int month1,int day1,int year2,int month2,int day2,bool *matherror); //计算两个日期之间的天数
int countdate(int *y,int *m,int *d,double n,bool *matherror); //计算某天之前/之后多少天的日期，返回该日期的星期数
inline bool checkyear(const int y,const int m,const int d); //判断某个日期是否合理
inline bool leepyear(const int year); //判断闰年
/*----------------------------------------------------------------------------------------------------------*/
double det(const int n,const double number[ ])
{
	double **a=new double*[n+2],b,temp;
	int i,j=0,k=1,m=1;
	bool t;
	for(i=0;i<=n+1;i++)
		a[i]=new double[n+2];
	for(i=1;i<=n*n;i++)
	{
		j++;
		if(j>n)
		{
			j=1;
			k++;
		}
		a[k][j]=number[i];
	}
	for(i=1;i<=n;i++)
	{
		t=(a[i][i]!=0);
		if(!t)
		{
			for(j=i+1;j<=n;j++)
				if(a[j][i]!=0)
				{
					t=1;
					for(k=1;k<=n;k++)
					{
						temp=a[i][k];
						a[i][k]=a[j][k];
						a[j][k]=temp;
					}
					m*=-1;
					break;
				}
			if(!t)
			{
				for(i=0;i<=n+1;i++)
		            delete []a[i];
	            delete []a;
				return 0;
			}
		}
		for(j=i+1;j<=n;j++)
		{
			b=-a[j][i]/a[i][i];
			for(k=i;k<=n;k++)
				a[j][k]+=b*a[i][k];
		}
	}
	b=m;
	for(i=1;i<=n;i++)
		b*=a[i][i];
	for(i=0;i<=n+1;i++)
		delete []a[i];
	delete []a;
	return b;
}
/*----------------------------------------------------------------------------------------------------------*/
double tr(const int n,const double number[ ])
{
	double r=0;
	int i;
	for(i=1;i<=n*n;i+=n+1)
		r+=number[i];
	return r;
}
/*----------------------------------------------------------------------------------------------------------*/
double rank(const int r,const int c,const double number[ ])
{
	double **a=new double*[r+2],b,temp;
	int i,j=0,k=1,l,rk=0;
	bool t;
	for(i=0;i<=r+1;i++)
		a[i]=new double[c+2];
    for(i=1;i<=r*c;i++)
	{
		j++;
		if(j>c)
		{
			j=1;
			k++;
		}
		a[k][j]=number[i];
	}
	t=1;
	j=0;
	for(i=1;i<=r;i++)
	{
		do
		{
			j++;
		    if(a[i][j]!=0)t=1;
		    else
			{
			    t=0;
			    for(k=i+1;k<=r;k++)
				    if(a[k][j]!=0)
					{
					    t=1;
					    for(l=j;l<=c;l++)
						{
						    temp=a[i][l];
						    a[i][l]=a[k][l];
						    a[k][l]=temp;
						}
					    break;
					}
			}
		}while(j<c&&!t);
		if(j>=c&&!t)break;
		for(k=i+1;k<=r;k++)
		{
			b=-a[k][j]/a[i][j];
			for(l=j;l<=c;l++)
				a[k][l]+=b*a[i][l];
		}
		if(j>=c)break;
	}
	t=0;
	for(i=r;i>=1;i--)
	{
		for(j=c;j>=1;j--)
		    if(fabs(a[i][j])>1e-14)
			{
			    rk=i;
				t=1;
			    break;
			}
		if(t)break;
	}
	for(i=0;i<r+1;i++)
		delete []a[i];
	delete []a;
	return minnum(rk,minnum(r,c));
}
/*----------------------------------------------------------------------------------------------------------*/
double fact(const double x,bool *matherror)
{
	int t=1;
	__int64 i,j;
	double r;
	if(x<0)t=-1;
	if(!isint(x))
	{
		*matherror=1;
		return 0;
	}
	i=__int64(fabs(x));
	if(i==0)return 1;
	r=1;
	for(j=2;j<=i;j++)
		r*=j;
	r*=t;
	return r;
}
/*----------------------------------------------------------------------------------------------------------*/
double factdouble(const double x,bool *matherror)
{
	int t=1;
	__int64 i,j;
	double r;
	if(x<0)t=-1;
	if(!isint(x))
	{
		*matherror=1;
		return 0;
	}
	i=__int64(fabs(x));
	if(i==0)return 1;
	r=1;
	for(j=i;j>=2;j-=2)
		r*=j;
	r*=t;
	return r;
}
/*----------------------------------------------------------------------------------------------------------*/
double count(double x,const double y,const int z,bool *matherror,const bool degreemode)
{
	if(*matherror==1)return 0;
	if(z==11||z==37||z==38||z==39||z==96||z==125)return x;
	if(z==1)return (x+y);
	if(z==2)return (x-y);
	if(z==3||z==103)return (x*y);
	if(z==4)
	{
		if(y==0)
		{
			*matherror=1;
			return 0;
		}
		return (x/y);
	}
	if(z==6)
	{
		if(fabs(y)<1e-15)
		{
			*matherror=1;
			return 0;
		}
		return __int64(x/y);
	}
	if(z==7)
	{
		if(fabs(y)<1e-15)
		{
			*matherror=1;
			return 0;
		}
		if(fabs(fabs(x)-fabs(__int64(x/y))*fabs(y))<1e-14)return 0;
		return fabs(x)-fabs(__int64(x/y))*fabs(y);
	}
	if(z==5)return __pow(x,y,matherror);
	if(z==12)
	{
		if(x<=0)
		{
			*matherror=1;
			return 0;
		}
		return log(x);
	}
	if(z==13)return fabs(x);
	if(z==15)
	{
		if(x<0)
		{
			*matherror=1;
			return 0;
		}
		return sqrt(x);
	}
	if(z==16)
	{
		if(degreemode)x=x*PI/180;
		if(fabs(sin(x))<5e-16)return 0;
		return sin(x);
	}
	if(z==17)
	{
		if(degreemode)x=x*PI/180;
		if(fabs(cos(x))<5e-16)return 0;
		return cos(x);
	}
	if(z==18)
	{
		if(degreemode)x=x*PI/180;
		if(fabs(tan(x))>1.633177e+016)*matherror=1;
		return tan(x);
	}
	if(z==19)
	{
		if(x<=0)
		{
			*matherror=1;
			return 0;
		}
		return log10(x);
	}
	if(z==20)return pow(Ee,x);
	if(z==21)
	{
		if(degreemode)x=x*PI/180;
		return count(1,tan(x),4,matherror,degreemode);
	}
	if(z==22)
	{
		if(degreemode)return asin(x)*180/PI;
		return asin(x);
	}
	if(z==23)
	{
		if(degreemode)return acos(x)*180/PI;
		return acos(x);
	}
	if(z==24)
	{
		if(degreemode)return atan(x)*180/PI;
		return atan(x);
	}
	if(z==25)return sinh(x);
	if(z==26)return cosh(x);
	if(z==27)return tanh(x);
	if(z==28)
	{
		if(x<=0||y<=0)
		{
			*matherror=1;
			return 0;
		}
		return log(y)/log(x);
	}
	if(z==29)return C(x,y,matherror);
	if(z==30)return P(x,y,matherror);
	if(z==31)return gcd(x,y,matherror);
	if(z==32)return lcm(x,y,matherror);
	if(z==34)return gauss(x);
	if(z==35)return x-gauss(x);
	if(z==36)return ran(x);
	if(z==45)return logiand(x,y);
	if(z==46)return logior(x,y);
	if(z==47)return logixor(x,y);
	if(z==55)
	{
		if(degreemode)x=x*PI/180;
		return 1/cos(x);
	}
	if(z==56)
	{
		if(degreemode)x=x*PI/180;
		return 1/sin(x);
	}
	if(z==57)return log(x+sqrt(x*x+1));
	if(z==58)return log(x+sqrt(x*x-1));
	if(z==59)return log((1+x)/(1-x))/2;
	if(z==71)return fib(x,matherror);
	if(z==72)return sgn(x);
	if(z==75)return int(x)&int(y);
	if(z==76)return int(x)|int(y);
	if(z==77)return int(x)^int(y);
	if(z==79)
	{
		if(!isint(y))*matherror=1;
		return int(x)<<int(y);
	}
	if(z==80)
	{
		if(!isint(y))*matherror=1;
		return int(x)>>int(y);
	}
	if(z==81)return round(x,y);
	if(z==82)
	{
		if(x>=0&&x<=4.3)return round(x,2);
		x=gauss(x);
	    if(x>100||x<0)
		{
		    *matherror=1;
		    return 0;
		}
		if(y==0)
		{
		    if(x<60)return 0;
		    return round(4-3*(100-x)*(100-x)/1600,2);
		}
		if(y==1)
		{
			if(x>=90)return 4;
			if(x>=80)return 3;
			if(x>=70)return 2;
			if(x>=60)return 1;
			return 0;
		}
		if(y==2)
		{
			if(x>=90)return 4;
			if(x>=85)return 3.7;
			if(x>=82)return 3.3;
			if(x>=78)return 3;
			if(x>=75)return 2.7;
			if(x>=72)return 2.3;
			if(x>=68)return 2;
			if(x>=64)return 1.5;
			if(x>=60)return 1;
			return 0;
		}
		if(y==3)return x/25;
		if(y==4)
		{
			if(x>=90)return 4.3;
			if(x>=85)return 4;
			if(x>=80)return 3.7;
			if(x>=75)return 3.3;
			if(x>=70)return 3;
			if(x>=65)return 2.7;
			if(x>=60)return 2;
			return 0;
		}
	}
	if(z==83)
	{
		if(degreemode)return (PI/2-atan(x))*180/PI;
		return PI/2-atan(x);
	}
	if(z==95)return 1-logixor(x,y);
	if(z==118)return gauss(x);
	if(z==119)return ceil(x);
	if(z==129)return factornumber(x);
	return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
inline bool cmp(const double x,const double y,const int z)
{
	if((z==51||z==52||z==53)&&fabs(x-y)<1e-15)return 1;
	if((z==49||z==52)&&x<y&&fabs(x-y)>1e-15)return 1;
	if((z==50||z==53)&&x>y&&fabs(x-y)>1e-15)return 1;
	return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
double gcd(double x,double y,bool *matherror)
{
	if(x<0)x*=-1;
	if(y<0)y*=-1;
	if(minnum(x,y)<1e-5)
	{
		*matherror=1;
		return 0;
	}
	double temp;
	if(x<y)
	{
		temp=x;
		x=y;
		y=temp;
	}
	temp=gauss(x/y)*y;
	if(fabs(x-temp)<1e-14)return y;
	return gcd(x-temp,y,matherror);
}
/*----------------------------------------------------------------------------------------------------------*/
int gcdn(int x,int y)
{
	if(x<y)
	{
		x^=y;
		y^=x;
		x^=y;
	}
	if(x%y==0)return y;
	return gcdn(x-(x/y)*y,y);
}
/*----------------------------------------------------------------------------------------------------------*/
double lcm(const double x,const double y,bool *matherror)
{
	return x*y/gcd(x,y,matherror);
}
/*----------------------------------------------------------------------------------------------------------*/
double ran(const double x)
{
	bool t=x<0?1:0;
	int i,b;
	double c=1,r2=0,a=1,r=0,m16[10]={1};
	for(i=1;i<=9;i++)
		m16[i]=16*m16[i-1];
	a=fabs(__int64(x));
	if(a>10000000000.0)a=10000000000.0;
	if(a<=1)
	{
		for(i=1;i<=8;i++)
		{
			c/=10;
			do
			{
				b=rand();
			}while(b>=32700);
			b%=10;
			r2+=b*c;
		}
		if(t)
		{
			b=rand()%2;
			if(b==1)r2*=-1;
		}
		return r2;
	}
	do
	{
		r=0;
	    for(i=0;i<9;i++)
		{
		    b=rand()%16;
		    r+=b*m16[i];
		}
	}while(r>=gauss(m16[9]/a)*a);
	r-=gauss(r/a)*a;
	if(t)
	{
		b=rand()%2;
		if(b==1)r*=-1;
	}
	return r;
}
/*----------------------------------------------------------------------------------------------------------*/
inline double minnum(const double x,const double y)
{
	if(x<y)return x;
	return y;
}
/*----------------------------------------------------------------------------------------------------------*/
inline double maxnum(const double x,const double y)
{
	if(x>y)return x;
	return y;
}
/*----------------------------------------------------------------------------------------------------------*/
double P(const double x,const double y,bool *matherror)
{
	if(x<0||y<0||!isint(x)||!isint(y))
	{
		*matherror=1;
		return 0;
	}
	__int64 i=__int64(x),j=__int64(y),tempint=i-j;
	double k=1;
	if(i<j)
	{
		*matherror=1;
		return 0;
	}
	for( ;i>tempint;i--)
		k*=i;
	return k;
}
/*----------------------------------------------------------------------------------------------------------*/
double C(const double x,const double y,bool *matherror)
{
	if(x<0||y<0||!isint(x)||!isint(y))
	{
		*matherror=1;
		return 0;
	}
	__int64 i=__int64(x),j=__int64(y);
	double k=1;
	if(i<j)
	{
		*matherror=1;
		return 0;
	}
	if(j>i-j)j=i-j;
	__int64 tempint=i-j;
	for( ;i>tempint;i--)
		k*=i;
	for(i=2;i<=j;i++)
		k/=i;
	return k;
}
/*----------------------------------------------------------------------------------------------------------*/
double gauss(const double x)
{
	double r=floor(x);
	if(fabs(r+1-x)<1e-12)r++;
	return r;
}
/*----------------------------------------------------------------------------------------------------------*/
inline bool isint(const double x)
{
	if(fabs(gauss(x)-x)<1e-12)return 1;
	return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
double __pow(const double x,double y,bool *matherror)
{
	if(x>=0)
	{
		if(x==0&&y<0)*matherror=1;
		return pow(x,y);
	}
	if(isint(y))return pow(x,gauss(y));
	int i,t=0,t2=0;
	if(y<0)
	{
		t=1;
		y*=-1;
	}
	double inty=gauss(y),tempy1,halftempy;
	y-=inty;
	for(i=2;i<=1000000;i++)
	{
		tempy1=y*i;
		if(fabs(tempy1-gauss(tempy1))<1e-10)
		{
			if(i%2==0)
			{
				*matherror=1;
				return 0;
			}
			else
			{
				halftempy=tempy1/2;
				if(fabs(halftempy-gauss(halftempy))<1e-10)t2=1;
				tempy1=pow(x,inty)*pow(-x,y);
				if(t==1)tempy1=1/tempy1;
				if(t2==0)tempy1*=-1;
				return tempy1;
			}
		}
	}
	*matherror=1;
	return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
inline double loginot(const double x)
{
	if(fabs(x)<1e-15)return 1;
	return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
inline double logiand(const double x,const double y)
{
	if((fabs(x)<1e-15)||(fabs(y)<1e-15))return 0;
	return 1;
}
/*----------------------------------------------------------------------------------------------------------*/
inline double logior(const double x,const double y)
{
	if((fabs(x)<1e-15)&&(fabs(y)<1e-15))return 0;
	return 1;
}
/*----------------------------------------------------------------------------------------------------------*/
inline double logixor(const double x,const double y)
{
	if((fabs(x)<1e-15)&&(fabs(y)<1e-15))return 0;
	if((fabs(x)<1e-15)||(fabs(y)<1e-15))return 1;
	return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
double fib(double x,bool *matherror)
{
	if(!isint(x))
	{
		*matherror=1;
		return 0;
	}
	x=gauss(x);
	double a1=(1+sqrt(5))/2,a2=(1-sqrt(5))/2;
	return (__pow(a1,x,matherror)-__pow(a2,x,matherror))/sqrt(5);
}
/*----------------------------------------------------------------------------------------------------------*/
inline short sgn(const double x)
{
	if(x==0)return 0;
	if(x>0)return 1;
	return -1;
}
/*----------------------------------------------------------------------------------------------------------*/
double round(double x,const double y)
{
	int bit=int(gauss(y)),t=1;
	if(bit<0||bit>307||fabs(x)>1e+18)return x;
	if(x<0)
	{
		t=-1;
		x=-x;
	}
	ofstream outfile;
	outfile.open("temp.txt",ios::out);
	outfile<<setiosflags(ios::fixed)<<setprecision(bit)<<x;
	outfile.close();
	ifstream infile;
	infile.open("temp.txt",ios::in);
	infile>>x;
	infile.close();
	remove("temp.txt");
	if(t==-1)x=-x;
	return x;
}
/*----------------------------------------------------------------------------------------------------------*/
double factornumber(const double x)
{
	if(!isint(x)||fabs(x)>=1e+15||fabs(x)<1e-10)return 0;
	__int64 n=__int64(gauss(fabs(x))),i,k=__int64(sqrt(n));
	if(n==1)return 1;
	int p;
	double r=1;
	for(i=2;i<=k;i++)
	{
		p=1;
		while(n%i==0)
		{
			p++;
			n/=i;
		}
		r*=p;
		if(n==1)break;
	}
	if(n!=1)r*=2;
	return r;
}
/*----------------------------------------------------------------------------------------------------------*/
double planeangle(const double x,const double y)
{
	if(y==0&&x>=0)return 0;
	if(x==0&&y>0)return PI/2;
	if(x==0&&y<0)return 3*PI/2;
	if(y==0&&x<0)return PI;
	int i;
	double a;
	if(x>0&&y>0)i=1;
	else if(x<0&&y>0)i=2;
	else if(x<0&&y<0)i=3;
	else i=4;
	a=atan(fabs(y)/fabs(x));
	if(i==1)return a;
	else if(i==2)return PI-a;
	else if(i==3)return PI+a;
	else if(i==4)return 2*PI-a;
	else return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
inline bool checkyear(const int y,const int m,const int d)
{
	if(y>=100000||y<=0||m<=0||m>=13||d<=0||d>=32)return 0;
	const int dayofmonth[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
	if(m!=2&&d>dayofmonth[m])return 0;
	if(m==2&&leepyear(y)&&d>29)return 0;
	if(m==2&&!leepyear(y)&&d>28)return 0;
	return 1;
}
/*----------------------------------------------------------------------------------------------------------*/
inline bool leepyear(const int year)
{
	if(year%400==0)return 1;
	if(year%4==0&&year%100!=0)return 1;
	return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
int daysbetween(int year1,int month1,int day1,int year2,int month2,int day2,bool *matherror)
{
	if(*matherror||!checkyear(year1,month1,day1)||!checkyear(year2,month2,day2))
	{
		*matherror=1;
		return 0;
	}
	bool t=1;
	if(year2<year1||(year2==year1&&month2<month1)||(year2==year1&&month2==month1&&day2<=day1))
	{
		year1^=year2;
		year2^=year1;
		year1^=year2;
		month1^=month2;
		month2^=month1;
		month1^=month2;
		day1^=day2;
		day2^=day1;
		day1^=day2;
		t=0;
	}
	if(year1==year2&&month1==month2)return t?day1-day2:day2-day1;
	const int dayofmonth[13]={31,31,28,31,30,31,30,31,31,30,31,30,31};
	int i,k=0;
	if(year1==year2)
	{
		for(i=month1+1;i<month2;i++)
			k+=dayofmonth[i];
		k+=dayofmonth[month1]-day1+day2;
		if(leepyear(year1)&&month1<=2&&month2>2)k++;
		return t?-k:k;
	}
	for(i=month1+1;i<=12;i++)
		k+=dayofmonth[i];
	k+=dayofmonth[month1]-day1;
	if(month1<=2&&leepyear(year1))k++;
	for(i=1;i<month2;i++)
		k+=dayofmonth[i];
	k+=day2;
	if(month2>2&&leepyear(year2))k++;
	k+=365*(year2-year1-1);
	if(year2-year1<=1000)
	{
		for(i=year1+1;i<year2;i++)
		    if(leepyear(i))k++;
	}
	else
	{
		int year1m4,year2m4,year1m100,year2m100,year1m400,year2m400;
    	for(year1m4=year1+1; ;year1m4++)
	    	if(year1m4%4==0)break;
	    for(year1m100=year1+1; ;year1m100++)
		    if(year1m100%100==0)break;
	    for(year1m400=year1+1; ;year1m400++)
		    if(year1m400%400==0)break;
	    for(year2m4=year2-1; ;year2m4--)
		    if(year2m4%4==0)break;
	    for(year2m100=year2-1; ;year2m100--)
		    if(year2m100%100==0)break;
	    for(year2m400=year2-1; ;year2m400--)
		    if(year2m400%400==0)break;
		k+=(year2m4-year1m4)/4-(year2m100-year1m100)/100+(year2m400-year1m400)/400+1;
	}
	return t?-k:k;
}
/*----------------------------------------------------------------------------------------------------------*/
int countdate(int *y,int *m,int *d,double pastday,bool *matherror)
{
	if(fabs(pastday)>36523883||!checkyear(*y,*m,*d)||!isint(pastday))
	{
		*matherror=1;
		return 0;
	}
	int i,temp,n=int(gauss(pastday));
	const int dayofmonth[13]={31,31,28,31,30,31,30,31,31,30,31,30,31};
	if(n==0)
	{
		i=daysbetween(*y,*m,*d,1,1,1,matherror);
		return (i%7+1)%7;
	}
	if(n>0)
	{
		if(n>1000)
		{
			temp=n/365+(*y);
			while(-daysbetween(*y,*m,*d,temp,1,1,matherror)>n)
				temp--;
			n+=daysbetween(*y,*m,*d,temp,1,1,matherror);
			*y=temp;
			*m=*d=1;
		}
		for(i=1;i<=n;i++)
		{
			(*d)++;
			if((*m)==2&&leepyear(*y))temp=1;
			else temp=0;
			if((*d)>(dayofmonth[(*m)]+temp))
			{
				*d=1;
				(*m)++;
			}
			if(*m==13)
			{
				*m=1;
				(*y)++;
			}
		}
	}
	else if(n<0)
	{
		n*=-1;
		if(n>1000)
		{
			temp=(*y)-n/365;
			while(daysbetween(*y,*m,*d,temp,1,1,matherror)>n)
				temp++;
			n-=daysbetween(*y,*m,*d,temp,1,1,matherror);
			*y=temp;
			*m=*d=1;
		}
		for(i=1;i<=n;i++)
		{
			(*d)--;
			if((*m)==3&&leepyear(*y))temp=1;
			else temp=0;
			if((*d)==0)
			{
				*d=dayofmonth[(*m)-1]+temp;
				(*m)--;
			}
			if((*m)==0)
			{
				*m=12;
				(*y)--;
			}
		}
	}
	if(!checkyear(*y,*m,*d))
	{
		*y=*m=*d=1;
		*matherror=1;
	}
	return 0;
}