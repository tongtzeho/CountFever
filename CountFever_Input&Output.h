using namespace std;
double inputvardata(double *(p[ ]),const double d,const bool t1); //输入变量
bool fracout(double x,const int form,const int blank,const int prec); //将某个数用分式或简单的有理根式输出（如果可能的话）
bool complexout(const double x,double y); //以分式或简单的有理根式输出复数
void complexdecout(const double x,const double y); //以小数输出复数
void diffout(string str,const bool *fracoutmode); //输出导函数
/*----------------------------------------------------------------------------------------------------------*/
double inputvardata(double *(p[ ]),const double d,const bool t1)
{
	char tempnum[500],tempnum2[500];
	int i,j,length;
	bool t=0;
	double a,b;
	do
	{
		cin.getline(tempnum,500,'\n');
		t=0;
		length=strlen(tempnum);
		if(t1&&length==0)return d;
		if(tempnum[0]=='M'&&tempnum[1]>='A'&&tempnum[1]<='Z') a=*p[tempnum[1]-64];
	    else if(tempnum[0]=='A'&&tempnum[1]=='N'&&tempnum[2]=='S') a=*p[0];
		else a=atof(tempnum);
		for(i=0;i<length;i++)
			if(tempnum[i]=='/')
			{
				t=1;
				break;
			}
		if(t)
		{
		    for(j=i+1;j<length;j++)
			    tempnum2[j-i-1]=tempnum[j];
		    if(tempnum2[0]=='M'&&tempnum2[1]>='A'&&tempnum2[1]<='Z') b=*p[tempnum2[1]-64];
	        else if(tempnum2[0]=='A'&&tempnum2[1]=='N'&&tempnum2[2]=='S') b=*p[0];
			else b=atof(tempnum2);
			if(b==0)
			{
				length=0;
				cout<<"出错！除数不能为零！重新输入 :\n";
			}
		}
	}while(length==0);
	if(t)return a/b;
	return a;
}
/*----------------------------------------------------------------------------------------------------------*/
bool fracout(double x,const int form,const int blank,const int prec)
{
	if(form==0)
	{
		cout<<setprecision(12)<<x;
		return 0;
	}
	int t=1,i,j;
	double numerator,tempx,eps=form==5?1e-7:1e-10;
	char temp[40];
	if(x<0)t=-1;
	x=fabs(x);
	if(x>10000||x<1e-4||fabs(x-gauss(x))<1e-10)
	{
		if(form==5)return 0;
		if(form==1)cout<<setw(blank)<<setprecision(prec)<<t*x;//1
		else if(form==2||form==3)
		{
			if(fabs(x)<1e-16)cout<<0;
			else cout<<setprecision(prec)<<t*x;
		}
		return 0;
	}
	ofstream outfile;
	ifstream infile;
	for(i=2;i<=20000;i++)
	{
		numerator=x*i;
		if(fabs(numerator-gauss(numerator))<1e-10)
		{
			if(form!=4)
			{
				outfile.open("temp.txt",ios::out);
				if(t==-1)outfile<<"-";
				outfile<<setprecision(9);
				if(form==5)outfile<<"(";
				outfile<<gauss(numerator)<<"/"<<i;
				if(form==5)outfile<<")";
				outfile<<endl;
				outfile.close();
				if(form==5)return 1;
				infile.open("temp.txt",ios::in);
				infile.getline(temp,40,'\n');
				infile.close();
				remove("temp.txt");
				if(form==1)
				{
					if(strlen(temp)<blank)cout<<setw(blank)<<temp;
					else cout<<setw(blank)<<setprecision(prec)<<t*x;
				}
				else cout<<temp;
			}
			else
			{
				if(t==-1)cout<<"-";
			    cout<<gauss(numerator)<<"/"<<i<<endl;
			}
			return 1;
		}
	}
	const short n[121]={2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,33,34,35,37,38,39,41,42,43,46,47,51,53,55,57,58,59,61,62,65,66,67,69,70,71,73,74,77,78,79,82,83,85,86,87,89,91,93,94,95,97,101,102,103,105,106,107,109,110,111,113,114,115,118,119,122,123,127,129,130,131,133,134,137,138,139,141,142,143,145,146,149,151,154,155,157,158,159,161,163,165,166,167,170,173,174,177,178,179,181,182,183,185,186,187,190,191,193,194,195,197,199};
	const double mathcon[7]={PI,Ee,log(2),PI*PI,1/PI,sqrt(PI),log(3)};
	const string mathconexp[7]={"π","e","ln(2)","π^2","π^-1","√π","ln(3)"};
	const string mathconexp_form5[7]={"pi","e","ln(2)","(pi^2)","(pi^-1)","sq(pi)","ln(3)"};
	for(i=0;i<=127;i++)
	{
		if(i<121)tempx=x/sqrt(n[i]);
		else tempx=x/mathcon[i-121];
		if(fabs(tempx-gauss(tempx))<1e-10)
		{
			if(form!=4)
			{
				outfile.open("temp.txt",ios::out);
				if(t==-1)outfile<<"-";
				outfile<<setprecision(9);
			    if(gauss(tempx)!=1)
				{
					if(form==5)outfile<<"(";
					outfile<<gauss(tempx);
				}
			    if(i<121)
				{
					if(form!=5)outfile<<"√"<<n[i]<<"\n";
					else if(gauss(tempx)==1)outfile<<"sq("<<n[i]<<")\n";
					else outfile<<"*sq("<<n[i]<<"))\n";
				}
			    else 
				{
					if(form!=5)outfile<<mathconexp[i-121]<<endl;
					else if(gauss(tempx)!=1)outfile<<"*"<<mathconexp_form5[i-121]<<")\n";
					else outfile<<mathconexp_form5[i-121]<<endl;
				}
				outfile.close();
				if(form==5)return 1;
				infile.open("temp.txt",ios::in);
				infile.getline(temp,40,'\n');
				infile.close();
				remove("temp.txt");
				if(form==1)
				{
					if(strlen(temp)<blank)cout<<setw(blank)<<temp;
					else cout<<setw(blank)<<setprecision(prec)<<t*x;
				}
				else cout<<temp;
			}
			else
			{
				if(t==-1)cout<<"-";
			    if(gauss(tempx)!=1)cout<<gauss(tempx);
			    if(i<121)cout<<"√"<<n[i]<<"\n";
			    else cout<<mathconexp[i-121]<<endl;
			}
			return 1;
		}
		for(j=2;j<=100;j++)
		{
			if(form==5)break;
			numerator=tempx*j;
			if(fabs(numerator-gauss(numerator))<1e-10)
			{
				if(form!=4)
				{
					outfile.open("temp.txt",ios::out);
					if(t==-1)outfile<<"-";
					outfile<<setprecision(9);
				    if(gauss(numerator)!=1)outfile<<gauss(numerator);
				    if(i<121)outfile<<"√"<<n[i]<<"/"<<j<<endl;
				    else outfile<<mathconexp[i-121]<<"/"<<j<<endl;
					outfile.close();
				    infile.open("temp.txt",ios::in);
				    infile.getline(temp,40,'\n');
				    infile.close();
					remove("temp.txt");
			    	if(form==1)
					{
					    if(strlen(temp)<blank)cout<<setw(blank)<<temp;
					    else cout<<setw(blank)<<setprecision(prec)<<t*x;
					}
					else cout<<temp;
				}
				else
				{
				    if(t==-1)cout<<"-";
				    if(gauss(numerator)!=1)cout<<gauss(numerator);
				    if(i<121)cout<<"√"<<n[i]<<"/"<<j<<endl;
				    else cout<<mathconexp[i-121]<<"/"<<j<<endl;
				}
				return 1;
			}
		}
	}
	if(form==5)return 0;
	if(form==3||form==4)
	{
		x*=t;
		t=1;
		int k,l,st,ed,g;
		double tempx2,tempx3;
		for(i=1;i<=40;i++)
		{
			if(i==1)
			{
				st=41;
				ed=60;
			}
			else if(i<=20)
			{
				st=21;
				ed=40;
			}
			else
			{
				st=1;
				ed=40;
			}
			for(t=1;t>=-1;t-=2)
				for(j=0;j<=17;j++)
					for(k=st;k<=ed;k++)
					{
					    if(j<=14)tempx=sqrt(n[j]);
					    else tempx=mathcon[j-15];
					    tempx2=x*i-t*tempx*k;
					    for(l=-1;l<=6;l++)
						{
						    if(l==-1)tempx3=tempx2;
						    else if(l<6)tempx3=tempx2/sqrt(n[l]);
						    else tempx3=tempx2/mathcon[l-6];
						    if(fabs(tempx3-gauss(tempx3))<1e-10)
							{
								g=gcdn(int(gauss(fabs(tempx3))),gcdn(i,k));
								if(form==3||i/g!=1)cout<<"(";
								if(int(gauss(fabs(tempx3)))/g!=1||l==-1)cout<<int(gauss(tempx3))/g;
								else if(int(gauss(tempx3))==-g)cout<<"-";
								if(l>=0&&l<6)cout<<"√"<<n[l];
								else if(l==6)cout<<mathconexp[l-6];
								if(t==1)cout<<"+";
								else cout<<"-";
								if(k/g!=1)cout<<k/g;
								if(j<15)cout<<"√"<<n[j];
								else cout<<mathconexp[j-15];
								if(form==3||i/g!=1)cout<<")";
								if(i/g!=1)cout<<"/"<<i/g;
								if(form==4)cout<<endl;
								return 1;
							}
						}
					}
		}
		t=1;
	}
	if(form==1)cout<<setw(blank)<<setprecision(prec)<<t*x;
	else if(form==2||form==3)cout<<setprecision(prec)<<t*x;
	return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
bool complexout(const double x,double y)
{
	if(fabs(x)<1e-14&&fabs(y)<1e-14)
	{
		cout<<"0";
		return 0;
	}
	bool t=y>0,froutx=0,frouty=0;
	if(fabs(x)<1e-14)
	{
		if(!t)
		{
			cout<<"-";
			y*=-1;
		}
		if(fabs(y-1)<1e-10)
		{
			frouty=0;
			cout<<"i";
		}
		else
		{
			frouty=fracout(y,3,NULL,12);
		    cout<<"i";
		}
		return frouty;
	}
	froutx=fracout(x,3,NULL,12);
	if(fabs(y)<1e-14)return froutx;
	if(t)cout<<"+";
	else
	{
		y*=-1;
		cout<<"-";
	}
	if(fabs(y-1)<1e-10)
	{
		frouty=0;
		cout<<"i";
	}
	else
	{
		frouty=fracout(y,3,NULL,12);
	    cout<<"i";
	}
	return froutx||frouty;
}
/*----------------------------------------------------------------------------------------------------------*/
void complexdecout(const double x,const double y)
{
	if(fabs(x)<1e-14)
	{
		if(fabs(y-1)<1e-10)cout<<"i";
		else if(fabs(y+1)<1e-10)cout<<"-i";
		else if(fabs(y)>1e-14)cout<<setprecision(12)<<y<<"i";
		else cout<<"0";
	}
	else 
	{
		cout<<setprecision(12)<<x;
		if(y>0&&fabs(y)>1e-14)cout<<"+";
		else if(fabs(y)>1e-14)cout<<"-";
		if(fabs(fabs(y)-1)<1e-10)cout<<"i";
		else if(fabs(y)>1e-14)cout<<setprecision(12)<<fabs(y)<<"i";
	}
}
/*----------------------------------------------------------------------------------------------------------*/
void diffout(string str,const bool *fracoutmode)
{
	int i,j,k,itemp,l,lstr;
	bool t;
	double tempnumber;
	string tempstr;
	if(*fracoutmode)
	{
		i=0;
		while(i<str.length())
		{
			if((str[i]>='0'&&str[i]<='9')||str[i]=='.')
			{
				itemp=i;
				k=0;
				while((str[itemp]>='0'&&str[itemp]<='9')||str[itemp]=='.')
				{
					itemp++;
					k++;
				}
				char tempch[30];
				for(j=0;j<30;j++)
					tempch[j]='\0';
				for(j=0;j<k;j++)
					tempch[j]=(str.substr(i,k))[j];
				tempnumber=atof(tempch);
				if(!fracout(tempnumber,5,NULL,NULL))
				{
					ofstream outfile;
	                outfile.open("temp.txt",ios::out);
	                outfile<<setprecision(10)<<tempnumber;
	                outfile.close();
				}
				ifstream infile;
                infile.open("temp.txt",ios::in);
	            infile>>tempstr;
	            infile.close();
	            remove("temp.txt");
				l=tempstr.length();
				lstr=str.length();
				str=str.substr(0,i)+tempstr+str.substr(itemp,lstr-itemp);
				i+=l;
			}
			else i++;
		}
	}
	int leftbr,rightbr,brsum=0,brleft,brtemp;
	bool *strdisp=new bool[str.length()+2],*brafind=new bool[str.length()+2];
	char leftch,rightch;
	for(i=0;i<=str.length();i++)
	{
		brafind[i]=0;
		strdisp[i]=1;
	}
	for(i=0;i<str.length();i++)
	{
		if(str[i]=='(')brsum++;
		if(str[i]=='('||str[i]==')')brafind[i]=1;
	}
	brleft=brsum;
	while(brleft>0)
	{
		for(i=0;i<str.length();i++)
		{
			if(str[i]=='('&&brafind[i])leftbr=i;
			else if(str[i]==')'&&brafind[i])
			{
				rightbr=i;
				break;
			}
		}
		brleft--;
		brafind[leftbr]=brafind[rightbr]=0;
		brtemp=t=0;
		for(i=leftbr+1;i<rightbr;i++)
		{
			if(str[i]=='('&&strdisp[i])brtemp++;
			if(str[i]==')'&&strdisp[i])brtemp--;
			if(brtemp==0&&(str[i]=='+'||str[i]=='-'||str[i]=='*'||str[i]=='/'||str[i]=='^'))
			{
				t=1;
				break;
			}
		}
		if(leftbr==0)leftch='(';
		else
		{
			for(j=leftbr-1;j>=0;j--)
				if(strdisp[j])
				{
					leftch=str[j];
					break;
				}
		}
		if(rightbr==str.length()-1)rightch=')';
		else
		{
			for(j=rightbr+1;j<str.length();j++)
				if(strdisp[j])
				{
					rightch=str[j];
					break;
				}
		}
		if(t)
		{
			if((leftch=='('||leftch=='+')&&(rightch==')'||rightch=='+'||rightch=='-'))strdisp[leftbr]=strdisp[rightbr]=0;
		}
		else
		{
			if(leftch=='('||leftch=='+'||leftch=='-'||leftch=='*'||leftch=='/'||leftch=='/')strdisp[leftbr]=strdisp[rightbr]=0;
		}
	}
	for(i=0;i<str.length();i++)
		if(strdisp[i])cout<<str[i];
	cout<<endl;
	delete []strdisp;
	delete []brafind;
}