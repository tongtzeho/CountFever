using namespace std;
struct countfeverdiffdata
{
	bool id,n[3];
	int cla,ope;
	double num[3],fac[3],nf[3];
	string f[3];
}diff[1000];
int diffsum,diffleft,leftb,rightb;
bool mathe=0,syntaxe=0;
void transf(char st[ ]); //翻译函数表达式
void differentiate(const int x,const int y,const int z); //在左右括号间求导(z为左括号对应的函数的代号)
string turndoubletostring(const double x); //将浮点型转成字符串
string turnfuntostring(const bool t,const double a,const string s,const double b,const double c,const int t0); //将函数表达式按一定法则转成字符串
bool issq(const string str); //判断是否根式
bool appequal(string st1,string st2); //判断两条表达式是否一样
int countdiff(const bool T1,const double A1,const string S1,const double B1,const double C1,const bool T2,const double A2,const string S2,const double B2,const double C2,bool *T,double *A,string *S,double *B,double *C,const int z); //对于各种运算进行求导
void diffdebugoutput(); //调试程序时的输出
/*----------------------------------------------------------------------------------------------------------*/
void transf(char st[ ])
{
	diffsum=diffleft=leftb=rightb=mathe=syntaxe=0;
	int i=0,j,length=strlen(st),t=0,t2=0,pointsum;
	char temp[500];
	while(i<length)
	{
		while(st[i]==' ')
			i++;
		if(st[i]=='+'||st[i]=='-')
		{
			if(t==0)
			{
				diffsum++;
				diff[diffsum].id=1;
				diff[diffsum].cla=4;
				if(st[i]=='+')diff[diffsum].ope=9;
				else diff[diffsum].ope=10;
			}
			if(t==1)
			{
				diffsum++;
				diff[diffsum].id=1;
				diff[diffsum].cla=2;
				if(st[i]=='+')diff[diffsum].ope=1;
				else diff[diffsum].ope=2;
				t=0;
			}
			t2=0;
		}
		else if(st[i]=='*'||st[i]=='/'||st[i]=='^')
		{
			diffsum++;
			diff[diffsum].id=1;
			diff[diffsum].cla=2;
			if(st[i]=='*')diff[diffsum].ope=3;
			else if(st[i]=='/')diff[diffsum].ope=4;
			else if(st[i]=='^')diff[diffsum].ope=5;
			t=t2=0;
		}
		else if(st[i]=='@'||st[i]=='#')
		{
			diffsum++;
			diff[diffsum].id=1;
			diff[diffsum].cla=3;
			if(st[i]=='@')diff[diffsum].ope=41;
			else if(st[i]=='#')diff[diffsum].ope=42;
			t=t2=1;
		}
		else if(st[i]=='s'||st[i]=='c'||st[i]=='('||st[i]=='l'||st[i]=='t'||st[i]=='a'||(st[i]=='e'&&st[i+1]=='x'))
		{
			if(t2==1)
			{
				diffsum++;
				diff[diffsum].id=1;
				diff[diffsum].cla=2;
				diff[diffsum].ope=3;
			}
			diffsum++;
			diff[diffsum].id=1;
			diff[diffsum].cla=5;
			if(st[i]=='s'&&st[i+1]=='q')
			{
				diff[diffsum].ope=15;
				if(st[i+2]=='r'&&st[i+3]=='t')i+=4;
				else i+=2;
			}
			else if(st[i]=='s'&&st[i+1]=='i'&&st[i+2]=='n')
			{
				diff[diffsum].ope=16;
				i+=3;
			}
			else if(st[i]=='l'&&st[i+1]=='n')
			{
				diff[diffsum].ope=12;
				i+=2;
			}
			else if(st[i]=='c'&&st[i+1]=='o'&&st[i+2]=='s')
			{
				diff[diffsum].ope=17;
				i+=3;
			}
			else if(st[i]=='l'&&st[i+1]=='g')
			{
				diff[diffsum].ope=19;
				i+=2;
			}
			else if(st[i]=='e'&&st[i+1]=='x')
			{
				diff[diffsum].ope=20;
				i+=3;
			}
			else if(st[i]=='t'&&(st[i+1]=='a'||st[i+1]=='g'))
			{
				diff[diffsum].ope=18;
				st[i+1]=='a'?i+=3:i+=2;
			}
			else if(st[i]=='c'&&(st[i+1]=='t'||(st[i+1]=='o'&&st[i+2]=='t')))
			{
				diff[diffsum].ope=21;
				i+=3;
			}
			else if(st[i]=='a'&&st[i+1]=='s'&&st[i+2]=='i')
			{
				diff[diffsum].ope=22;
				i+=4;
			}
			else if(st[i]=='a'&&st[i+1]=='c'&&st[i+2]=='o'&&st[i+3]=='s')
			{
				diff[diffsum].ope=23;
				i+=4;
			}
			else if(st[i]=='a'&&st[i+1]=='t'&&(st[i+2]=='a'||st[i+2]=='g'))
			{
				diff[diffsum].ope=24;
				st[i+2]=='a'?i+=4:i+=3;
			}
			else if(st[i]=='s'&&st[i+1]=='h')
			{
				diff[diffsum].ope=25;
				i+=2;
			}
			else if(st[i]=='c'&&st[i+1]=='h')
			{
				diff[diffsum].ope=26;
				i+=2;
			}
			else if(st[i]=='t'&&st[i+1]=='h')
			{
				diff[diffsum].ope=27;
				i+=2;
			}
			else if(st[i]=='s'&&st[i+1]=='e')
			{
				diff[diffsum].ope=55;
				i+=3;
			}
			else if(st[i]=='c'&&st[i+1]=='s')
			{
				diff[diffsum].ope=56;
				i+=3;
			}
			else if(st[i]=='a'&&st[i+1]=='s'&&st[i+2]=='h')
			{
				diff[diffsum].ope=57;
				i+=3;
			}
			else if(st[i]=='a'&&st[i+1]=='c'&&st[i+2]=='h')
			{
				diff[diffsum].ope=58;
				i+=3;
			}
			else if(st[i]=='a'&&st[i+1]=='t'&&st[i+2]=='h')
			{
				diff[diffsum].ope=59;
				i+=3;
			}
			else if(st[i]=='a'&&st[i+1]=='c'&&(st[i+2]=='t'||st[i+2]=='o'))
			{
				diff[diffsum].ope=83;
				i+=4;
			}
			else if(st[i]=='(')diff[diffsum].ope=11;
			else
			{
				syntaxe=1;
				return;
			}
			t=t2=0;
		}
		else if(st[i]==')')
		{
			diffsum++;
			diff[diffsum].id=1;
			diff[diffsum].cla=6;
			diff[diffsum].ope=14;
			t=t2=1;
		}
		else if((st[i]>='0'&&st[i]<='9')||st[i]=='.')
		{
			int templength=strlen(temp);
			pointsum=0;
			for(j=0;j<=templength;j++)
				temp[j]='\0';
			j=i;
			while((st[j]>=48&&st[j]<=57)||st[j]=='.')
			{
				temp[j-i]=st[j];
				if(st[j]=='.')pointsum++;
				if(pointsum>=2)
				{
					syntaxe=1;
					return;
				}
				j++;
			}
			i=j-1;
			diffsum++;
			diff[diffsum].id=1;
			diff[diffsum].cla=1;
			diff[diffsum].num[1]=atof(temp);
			diff[diffsum].num[2]=0;
			diff[diffsum].n[1]=diff[diffsum].n[2]=0;
			t=t2=1;
		}
		else if(st[i]=='e'||(st[i]=='p'&&st[i+1]=='i'))
		{
			if(t2==1)
			{
				diffsum++;
				diff[diffsum].id=1;
				diff[diffsum].cla=2;
				diff[diffsum].ope=3;
			}
			diffsum++;
			diff[diffsum].id=1;
			diff[diffsum].cla=1;
			diff[diffsum].num[2]=diff[diffsum].n[1]=diff[diffsum].n[2]=0;
			if(st[i]=='e')diff[diffsum].num[1]=Ee;
			else if(st[i]=='p'&&st[i+1]=='i')
			{
				diff[diffsum].num[1]=PI;
				i++;
			}
			t=t2=1;
		}
		else if(st[i]=='x'||st[i]=='X')
		{
			if(t2==1)
			{
				diffsum++;
				diff[diffsum].id=1;
				diff[diffsum].cla=2;
				diff[diffsum].ope=3;
			}
			diffsum++;
			diff[diffsum].id=1;
			diff[diffsum].cla=1;
			diff[diffsum].num[1]=0;
			diff[diffsum].num[2]=1;
			diff[diffsum].n[1]=1;
			diff[diffsum].n[2]=0;
			diff[diffsum].fac[1]=diff[diffsum].nf[1]=1;
			diff[diffsum].f[1]="X";
			t=t2=1;
		}
		else
		{
			syntaxe=1;
			return;
		}
		i++;
	}
	diffleft=diffsum;
	leftb=rightb=0;
	for(i=1;i<=diffsum;i++)
	{
		if(diff[i].cla==5)leftb++;
		else if(diff[i].cla==6)rightb++;
	}
	if(leftb>rightb)
	{
		for(i=diffleft+1;i<=diffleft+leftb-rightb;i++)
		{
			diffsum++;
			diff[diffsum].id=1;
			diff[diffsum].cla=6;
			diff[diffsum].ope=14;
		}
	}
	else if(leftb<rightb)
	{
		for(i=diffleft-1;i>diffleft-rightb+leftb;i--)
		{
			diffsum--;
			diff[diffsum].id=0;
		}
	}
	diffleft=diffsum;
}
/*----------------------------------------------------------------------------------------------------------*/
string turndoubletostring(const double x)
{
	string str;
	ofstream outfile;
	outfile.open("temp.txt",ios::out);
	outfile<<setprecision(14)<<x;
	outfile.close();
	ifstream infile;
	infile.open("temp.txt",ios::in);
	infile>>str;
	infile.close();
	remove("temp.txt");
	return str;
}
/*----------------------------------------------------------------------------------------------------------*/
string turnfuntostring(const bool t,const double a,const string s,const double b,const double c,const int t0)
{
	bool frat=s[0]=='/';
	if(t==0)
	{
		if(t0==2)return "("+turndoubletostring(c)+")";
		return turndoubletostring(c);
	}
	if(fabs(c)<1e-8)
	{
		if(frat)
		{
			if(fabs(b-1)<1e-8&&t0>=1)return "("+turndoubletostring(a)+s+")";
			else if(fabs(b-1)<1e-8)return turndoubletostring(a)+s;
			else
			{
				if(fabs(a-1)<1e-8)
				{
					if(t0>=1)return "(1"+s+"^"+turndoubletostring(b)+")";
					else return "1"+s+"^"+turndoubletostring(b);
				}
				else if(fabs(a+1)<1e-8)
				{
					if(t0>=1)return "(-1"+s+"^"+turndoubletostring(b)+")";
					else return "-1"+s+"^"+turndoubletostring(b);
				}
				else
				{
					if(t0>=1)return "("+turndoubletostring(a)+"*1"+s+"^"+turndoubletostring(b)+")";
					else return turndoubletostring(a)+"*1"+s+"^"+turndoubletostring(b);
				}
			}
		}
		if(fabs(a-1)<1e-8)
		{
			if(fabs(b-1)<1e-8&&t0==3&&s.length()==1)return s;
			else if(fabs(b-1)<1e-8&&t0>=1)return "("+s+")";
			else if(fabs(b-1)<1e-8)return s;
			else if(t0>=1)return "(("+s+")^"+turndoubletostring(b)+")";
			else return "("+s+")^"+turndoubletostring(b);
		}
		else if(fabs(a+1)<1e-8)
		{
			if(fabs(b-1)<1e-8&&t0>=1)return "(-"+s+")";
			else if(fabs(b-1)<1e-8)return "-"+s;
			else if(t0>=1)return "(-(("+s+")^"+turndoubletostring(b)+"))";
			else return "-(("+s+")^"+turndoubletostring(b)+")";
		}
		else
		{
			if(fabs(b-1)<1e-8&&t0>=1)return "("+turndoubletostring(a)+"*"+s+")";
			else if(fabs(b-1)<1e-8)return turndoubletostring(a)+"*"+s;
			else if(t0>=1)return "("+turndoubletostring(a)+"*("+s+")^"+turndoubletostring(b)+")";
			else return turndoubletostring(a)+"*("+s+")^"+turndoubletostring(b);
		}
	}
	else
	{
		string tempst=c>0?"+":"-";
		if(frat)
		{
			if(fabs(b-1)<1e-8)return "("+turndoubletostring(a)+s+tempst+turndoubletostring(fabs(c))+")";
			else
			{
				if(fabs(a-1)<1e-8)return "((1"+s+"^"+turndoubletostring(b)+")"+tempst+turndoubletostring(fabs(c))+")";
				else if(fabs(a+1)<1e-8)return "((-1"+s+"^"+turndoubletostring(b)+")"+tempst+turndoubletostring(fabs(c))+")";
				else return "("+turndoubletostring(a)+"*1"+s+"^"+turndoubletostring(b)+tempst+turndoubletostring(fabs(c))+")";
			}
		}
		if(fabs(a-1)<1e-8)
		{
			if(fabs(b-1)<1e-8)return "("+s+tempst+turndoubletostring(fabs(c))+")";
			else return "(("+s+")^"+turndoubletostring(b)+tempst+turndoubletostring(fabs(c))+")";
		}
		else if(fabs(a+1)<1e-8)
		{
			if(fabs(b-1)<1e-8)return "(-"+s+tempst+turndoubletostring(fabs(c))+")";
			else return "(-(("+s+")^"+turndoubletostring(b)+")"+tempst+turndoubletostring(fabs(c))+")";
		}
		else
		{
			if(fabs(b-1)<1e-8)return "("+turndoubletostring(a)+"*"+s+tempst+turndoubletostring(fabs(c))+")";
			else return "("+turndoubletostring(a)+"*("+s+")^"+turndoubletostring(b)+tempst+turndoubletostring(fabs(c))+")";
		}
	}
}
/*----------------------------------------------------------------------------------------------------------*/
void differentiate(const int x,const int y,const int z)
{
	int i,j,k,temppos1,temppos2;
	bool find;
	for(i=x+1;i<y;i++)
		if(diff[i].id==1&&diff[i].cla==3)
		{
			diff[i].id=0;
			diffleft--;
			find=0;
			for(j=i-1;j>x;j--)
				if(diff[j].id==1&&diff[j].cla==1)
				{
					find=1;
					break;
				}
			if(find==0)
			{
				syntaxe=1;
				return;
			}
			if(diff[j].cla==1)
			{
				k=diff[i].ope-39;
				string tempst[3];
			    double tempa[3],tempb[3],tempc[3];
			    bool tempt[3];
				temppos1=j;
				countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,k-1,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],5);
				countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
				countdiff(tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],0,0,"0",0,k,&diff[temppos1].n[2],&diff[temppos1].fac[2],&diff[temppos1].f[2],&diff[temppos1].nf[2],&diff[temppos1].num[2],3);
				countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,k,&diff[temppos1].n[1],&diff[temppos1].fac[1],&diff[temppos1].f[1],&diff[temppos1].nf[1],&diff[temppos1].num[1],5);
			}
		}
	for(i=x+1;i<y;i++)
		if(diff[i].id==1&&diff[i].cla==4)
		{
			diff[i].id=0;
			diffleft--;
			find=0;
			if(diff[i].ope==10)
			{
				for(j=i+1;j<y;j++)
					if(diff[j].id==1&&(diff[j].cla==1||diff[j].cla==4))
					{
						find=1;
						if(diff[j].cla==4)
						{
							if(diff[j].ope==9)diff[j].ope=10;
							else diff[j].ope=9;
						}
						else
						{
							diff[j].num[1]*=-1;
							diff[j].num[2]*=-1;
							for(k=1;k<=2;k++)
								if(diff[j].n[k])
									diff[j].fac[k]*=-1;
						}
						break;
					}
				if(!find)
				{
					syntaxe=1;
					return;
				}
			}
		}	
	for(i=x+1;i<y;i++)
		if(diff[i].id==1&&diff[i].cla==2&&diff[i].ope==5)
		{
			find=0;
			diff[i].id=0;
			for(j=i-1;j>x;j--)
				if(diff[j].id==1&&diff[j].cla==1)
				{
					find=1;
					temppos1=j;
					break;
				}
			if(!find)
			{
				syntaxe=1;
				return;
			}
			find=0;
			for(j=i+1;j<y;j++)
				if(diff[j].id==1&&diff[j].cla==1)
				{
					find=1;
					temppos2=j;
					diff[j].id=0;
					break;
				}
			if(!find)
			{
				syntaxe=1;
				return;
			}
			diffleft-=2;
			string tempst[4];
			double tempa[4],tempb[4],tempc[4];
			bool tempt[4];
			k=countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],diff[temppos2].n[1],diff[temppos2].fac[1],diff[temppos2].f[1],diff[temppos2].nf[1],diff[temppos2].num[1],&tempt[1],&tempa[1],&tempst[1],&tempb[1],&tempc[1],diff[i].ope);
			if(k==1||k==2||k==6||k==7)
			{
				tempa[2]=tempc[2]=tempt[2]=0;
				tempst[2]="0";
				tempb[2]=1;
			}
			else if(k==3)
			{
				tempt[2]=diff[temppos1].n[2];
				tempa[2]=diff[temppos1].fac[2];
				tempb[2]=diff[temppos1].nf[2];
				tempc[2]=diff[temppos1].num[2];
				tempst[2]=diff[temppos1].f[2];
			}
			else if(k==4||k==5||k==10)
			{
				countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,diff[temppos2].num[1]-1,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],5);
				countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
				countdiff(tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],0,0,"0",0,diff[temppos2].num[1],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
			}
			else if(k==8)
			{
				countdiff(diff[temppos2].n[2],diff[temppos2].fac[2],diff[temppos2].f[2],diff[temppos2].nf[2],diff[temppos2].num[2],0,0,"0",1,log(diff[temppos1].num[1]),&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],3);
				countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],tempt[1],tempa[1],tempst[1],tempb[1],tempc[1],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
			}
			else
			{
				countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],12);
				countdiff(diff[temppos2].n[2],diff[temppos2].fac[2],diff[temppos2].f[2],diff[temppos2].nf[2],diff[temppos2].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],3);
				countdiff(diff[temppos2].n[1],diff[temppos2].fac[1],diff[temppos2].f[1],diff[temppos2].nf[1],diff[temppos2].num[1],diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],4);
				countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
				countdiff(tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],1);
				countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],tempt[1],tempa[1],tempst[1],tempb[1],tempc[1],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
			}
			for(j=1;j<=2;j++)
			{
				diff[temppos1].n[j]=tempt[j];
				diff[temppos1].fac[j]=tempa[j];
				diff[temppos1].f[j]=tempst[j];
				diff[temppos1].nf[j]=tempb[j];
				diff[temppos1].num[j]=tempc[j];
			}
		}

	for(i=x+1;i<y;i++)
		if(diff[i].id==1&&diff[i].cla==2&&(diff[i].ope==3||diff[i].ope==4))
		{
			find=0;
			diff[i].id=0;
			for(j=i-1;j>x;j--)
				if(diff[j].id==1&&diff[j].cla==1)
				{
					find=1;
					temppos1=j;
					break;
				}
			if(!find)
			{
				syntaxe=1;
				return;
			}
			find=0;
			for(j=i+1;j<y;j++)
				if(diff[j].id==1&&diff[j].cla==1)
				{
					find=1;
					temppos2=j;
					diff[j].id=0;
					break;
				}
			if(!find)
			{
				syntaxe=1;
				return;
			}
			diffleft-=2;
			string tempst[4];
			double tempa[4],tempb[4],tempc[4];
			bool tempt[4];
			k=countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],diff[temppos2].n[1],diff[temppos2].fac[1],diff[temppos2].f[1],diff[temppos2].nf[1],diff[temppos2].num[1],&tempt[1],&tempa[1],&tempst[1],&tempb[1],&tempc[1],diff[i].ope);
			if(diff[i].ope==3)
			{
				if(k==1||k==2||k==4||k==6)
				{
					tempa[2]=tempc[2]=tempt[2]=0;
					tempst[2]="0";
					tempb[2]=1;
				}
				else if(k==3)countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],0,0,"0",1,diff[temppos2].num[1],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],diff[i].ope);
				else if(k==5)countdiff(diff[temppos2].n[2],diff[temppos2].fac[2],diff[temppos2].f[2],diff[temppos2].nf[2],diff[temppos2].num[2],0,0,"0",1,diff[temppos1].num[1],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],diff[i].ope);
				else if(k==7||k==9)
				{
					countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],diff[temppos2].n[1],diff[temppos2].fac[1],diff[temppos2].f[1],diff[temppos2].nf[1],diff[temppos2].num[1],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],3);
					countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],0,0,"0",1,2,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],diff[i].ope);
                }
				else
				{
					countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],diff[temppos2].n[1],diff[temppos2].fac[1],diff[temppos2].f[1],diff[temppos2].nf[1],diff[temppos2].num[1],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],3);
					countdiff(diff[temppos2].n[2],diff[temppos2].fac[2],diff[temppos2].f[2],diff[temppos2].nf[2],diff[temppos2].num[2],diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],3);
					countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],1);
				}
			}
			else
			{
				if(k==1||k==3||k==6||k==8)
				{
					tempa[2]=tempc[2]=tempt[2]=0;
					tempst[2]="0";
					tempb[2]=1;
				}
				else if(k==2)countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],0,0,"0",1,diff[temppos2].num[1],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],diff[i].ope);
				else if(k==4)
				{
					countdiff(diff[temppos2].n[2],diff[temppos2].fac[2],diff[temppos2].f[2],diff[temppos2].nf[2],diff[temppos2].num[2],0,0,"0",1,-diff[temppos1].num[1],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],3);
					countdiff(diff[temppos2].n[1],diff[temppos2].fac[1],diff[temppos2].f[1],diff[temppos2].nf[1],diff[temppos2].num[1],diff[temppos2].n[1],diff[temppos2].fac[1],diff[temppos2].f[1],diff[temppos2].nf[1],diff[temppos2].num[1],&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],3);
					countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
				}
				else
				{
					countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],diff[temppos2].n[1],diff[temppos2].fac[1],diff[temppos2].f[1],diff[temppos2].nf[1],diff[temppos2].num[1],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
					countdiff(diff[temppos2].n[2],diff[temppos2].fac[2],diff[temppos2].f[2],diff[temppos2].nf[2],diff[temppos2].num[2],diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],3);
					countdiff(tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],2);
					countdiff(diff[temppos2].n[1],diff[temppos2].fac[1],diff[temppos2].f[1],diff[temppos2].nf[1],diff[temppos2].num[1],0,0,"0",1,2,&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],5);
					countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
				}
			}
			for(j=1;j<=2;j++)
			{
				diff[temppos1].n[j]=tempt[j];
				diff[temppos1].fac[j]=tempa[j];
				diff[temppos1].f[j]=tempst[j];
				diff[temppos1].nf[j]=tempb[j];
				diff[temppos1].num[j]=tempc[j];
			}
		}
	for(i=x+1;i<y;i++)
		if(diff[i].id==1&&diff[i].cla==2&&(diff[i].ope==1||diff[i].ope==2))
		{
			find=0;
			diff[i].id=0;
			for(j=i-1;j>x;j--)
				if(diff[j].id==1&&diff[j].cla==1)
				{
					find=1;
					temppos1=j;
					break;
				}
			if(!find)
			{
				syntaxe=1;
				return;
			}
			find=0;
			for(j=i+1;j<y;j++)
				if(diff[j].id==1&&diff[j].cla==1)
				{
					find=1;
					temppos2=j;
					diff[j].id=0;
					break;
				}
			if(!find)
			{
				syntaxe=1;
				return;
			}
			diffleft-=2;
			for(j=1;j<=2;j++)
			{
				string tempst;
				double tempa,tempb,tempc;
				bool tempt;
				k=countdiff(diff[temppos1].n[j],diff[temppos1].fac[j],diff[temppos1].f[j],diff[temppos1].nf[j],diff[temppos1].num[j],diff[temppos2].n[j],diff[temppos2].fac[j],diff[temppos2].f[j],diff[temppos2].nf[j],diff[temppos2].num[j],&tempt,&tempa,&tempst,&tempb,&tempc,diff[i].ope);
				diff[temppos1].n[j]=tempt;
				diff[temppos1].fac[j]=tempa;
				diff[temppos1].f[j]=tempst;
				diff[temppos1].nf[j]=tempb;
				diff[temppos1].num[j]=tempc;
			}
		}
	if(z!=11)
	{
		temppos1=0;
		for(i=x+1;i<y;i++)
			if(diff[i].id==1&&diff[i].cla==1)
			{
				temppos1=i;
				break;
			}
		if(temppos1==0)
		{
			syntaxe=1;
			return;
		}
		string tempst[4];
		double tempa[4],tempb[4],tempc[4];
		bool tempt[4];
		k=countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[1],&tempa[1],&tempst[1],&tempb[1],&tempc[1],z);
        if(z==12)countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
		else if(z==15)
		{
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],0,0,"0",1,0.5,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],3);
			countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],tempt[1],tempa[1],tempst[1],tempb[1],tempc[1],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
		}
		else if(z==16)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],17);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
		}
		else if(z==17)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],16);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],3);
			countdiff(tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],0,0,"0",1,-1,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
		}
		else if(z==18)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],17);
			countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],0,0,"0",1,2,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],5);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
		}
		else if(z==19)
		{
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],4);
			countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],0,0,"0",1,log(10),&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
		}
		else if(z==20)countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[1],tempa[1],tempst[1],tempb[1],tempc[1],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
		else if(z==21)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],16);
			countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],0,0,"0",1,2,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],5);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],4);
			countdiff(tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],0,0,"0",1,-1,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
		}
		else if(z==22||z==23)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,2,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],5);
			countdiff(0,0,"0",0,1,tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],2);
			countdiff(tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],0,0,"0",1,0,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],15);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
			if(z==23)countdiff(tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],0,0,"0",1,-1,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
		}
		else if(z==24||z==83)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,2,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],5);
			countdiff(0,0,"0",0,1,tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],1);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
			if(z==83)countdiff(tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],0,0,"0",1,-1,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);		
		}
		else if(z==25||z==26)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],51-z);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
		}
		else if(z==27)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],26);
			countdiff(tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],0,0,"0",1,2,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],5);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
		}
		else if(z==55)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],16);
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],17);
			countdiff(tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],0,0,"0",1,2,&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],5);
			countdiff(tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],4);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
		}
		else if(z==56)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],17);
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,0,&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],16);
			countdiff(tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],0,0,"0",1,2,&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],5);
			countdiff(tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],4);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
			countdiff(tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],0,0,"0",1,-1,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],3);
		}
		else if(z==57||z==58)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,2,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],5);
			countdiff(tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],0,0,"0",0,1,&tempt[3],&tempa[3],&tempst[3],&tempb[3],&tempc[3],z-56);
			countdiff(tempt[3],tempa[3],tempst[3],tempb[3],tempc[3],0,0,"0",1,0,&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],15);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
		}
		else if(z==59)
		{
			countdiff(diff[temppos1].n[1],diff[temppos1].fac[1],diff[temppos1].f[1],diff[temppos1].nf[1],diff[temppos1].num[1],0,0,"0",1,2,&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],5);
			countdiff(0,0,"0",0,1,tempt[2],tempa[2],tempst[2],tempb[2],tempc[2],&tempt[0],&tempa[0],&tempst[0],&tempb[0],&tempc[0],2);
			countdiff(diff[temppos1].n[2],diff[temppos1].fac[2],diff[temppos1].f[2],diff[temppos1].nf[2],diff[temppos1].num[2],tempt[0],tempa[0],tempst[0],tempb[0],tempc[0],&tempt[2],&tempa[2],&tempst[2],&tempb[2],&tempc[2],4);
		}
		for(j=1;j<=2;j++)
		{
			diff[temppos1].n[j]=tempt[j];
			diff[temppos1].fac[j]=tempa[j];
			diff[temppos1].f[j]=tempst[j];
			diff[temppos1].nf[j]=tempb[j];
			diff[temppos1].num[j]=tempc[j];
		}
	}
}
/*----------------------------------------------------------------------------------------------------------*/
bool issq(const string str)
{
	int l=str.length(),i,ls=0,rs=0;
	if(!(str[0]=='s'&&str[1]=='q'&&str[2]=='r'&&str[3]=='t'))return 0;
	for(i=0;i<l;i++)
	{
		if(str[i]=='(')ls++;
		if(str[i]==')')rs++;
		if(ls==rs&&ls>0&&i<l-1)return 0;
	}
	return 1;
}
/*----------------------------------------------------------------------------------------------------------*/
bool appequal(string st1,string st2)
{
	if(st1==st2)return 1;
	int l1=st1.length(),l2=st2.length();
	if(l1==l2||(l1-l2)%2!=0)return 0;
	if(l1>l2)
	{
		int i;
		for(i=0;i<(l1-l2)/2;i++)
			if(st1[i]!='(')return 0;
		st1=st1.substr((l1-l2)/2,l2);
		return st1==st2;
	}
	else
	{
		int i;
		for(i=0;i<(l2-l1)/2;i++)
			if(st2[i]!='(')return 0;
		st2=st2.substr((l2-l1)/2,l1);
		return st1==st2;
	}
}
/*----------------------------------------------------------------------------------------------------------*/
int countdiff(const bool T1,const double A1,const string S1,const double B1,const double C1,const bool T2,const double A2,const string S2,const double B2,const double C2,bool *T,double *A,string *S,double *B,double *C,const int z)
{
	if(z==1||z==2)
	{
		*C=count(C1,C2,z,&mathe,NULL);
		if(T1&&T2)
		{
			if((S1==S2||appequal(S1,S2))&&fabs(B1-B2)<1e-8)
			{
				*A=count(A1,A2,z,&mathe,NULL);
				if(fabs(*A)<1e-8)
				{
					*T=0;
					*S="0";
					*A=0;
					*B=1;
					return 1;
				}
				else
				{
					*T=1;
					*S=S1;
					*B=B1;
					return 2;
				}
			}
			else
			{
				*T=1;
				*A=1;
				*B=1;
				*S="("+turnfuntostring(T1,A1,S1,B1,0,0);
				string tempst=turnfuntostring(T2,fabs(A2),S2,B2,0,0);
				if(tempst[0]=='-')
				{
					if(z==2)tempst[0]='+';
					*S+=tempst+")";
					return 3;
				}
				if((z==1&&A2>0)||(z==2&&A2<0))*S+="+"+tempst+")";
				else *S+="-"+tempst+")";
				return 3;
			}
		}
		else if(T1||T2)
		{
			*T=1;
			if(T1)
			{
				*A=A1;
				*B=B1;
				*S=S1;
				return 4;
			}
			else
			{
				*A=z==1?A2:-A2;
				*B=B2;
				*S=S2;
				return 5;
			}
		}
		else
		{
			*A=*T=0;
			*B=1;
			*S="0";
			return 6;
		}
	}
	if(z==3)
	{
		if(!T1&&!T2)
		{
			*A=*T=0;
			*S="0";
			*B=1;
			*C=C1*C2;
			return 1;
		}
		else if(T1+T2==1)
		{
			if(T1==1)
			{
				if(fabs(C2)<1e-8)
				{
					*A=*T=0;
					*S="0";
					*B=1;
					*C=0;
					return 2;
				}
				else
				{
					*T=1;
					*S=S1;
					*B=B1;
					*A=A1*C2;
					*C=C1*C2;
					return 3;
				}
			}
			else
			{
				if(fabs(C1)<1e-8)
				{
					*A=*T=0;
					*S="0";
					*B=1;
					*C=0;
					return 4;
				}
				else
				{
					*T=1;
					*S=S2;
					*B=B2;
					*A=A2*C1;
					*C=C1*C2;
					return 5;
				}
			}
		}
		else
		{
			if((S1==S2||appequal(S1,S2))&&fabs(B1-B2)<1e-8&&fabs(A1-A2)<1e-8&&fabs(C1-C2)<1e-8)
			{
				countdiff(1,A1,S1,B1,C1,0,0,"0",1,2,T,A,S,B,C,5);
				return 9;
			}
			else if(fabs(C1)<1e-8&&fabs(C2)<1e-8&&(S1==S2||appequal(S1,S2)))
			{
				if(fabs(B1+B2)<1e-8)
				{
					*A=*T=0;
					*S="0";
					*B=1;
					*C=A1*A2;
					return 6;
				}
				else
				{
					*T=1;
					*A=A1*A2;
					*B=B1+B2;
					*S=S1;
					*C=0;
					if(fabs(B1-B2)<1e-8)return 7;
					else return 8;
				}
			}
			else if(S2[0]=='/'&&fabs(C2)<1e-8)
			{
				int k,l=S2.length();
				k=countdiff(1,A1,S1,B1,C1,1,A2,S2.substr(1,l-1),B2,0,T,A,S,B,C,4);
				if(k==1||k==3||k==6||k==8)return 6;
				else return 10;
			}
			else if(S1[0]=='/'&&fabs(C1)<1e-8&&fabs(A1)!=0)
			{
				int k,l=S1.length();
				k=countdiff(1,A2,S2,B2,C2,1,1/A1,S1.substr(1,l-1),B1,0,T,A,S,B,C,4);
				if(k==1||k==3||k==6||k==8)return 6;
				else return 10;
			}
			else
			{
				*A=*B=*T=1;
				*C=0;
				*S=turnfuntostring(1,A1,S1,B1,C1,0)+"*"+turnfuntostring(1,A2,S2,B2,C2,0);
				return 10;
			}
		}
	}
	if(z==4)
	{
		if(!T1&&!T2)
		{
			*A=*T=0;
			*S="0";
			*B=1;
			if(C2==0)mathe=1;
			*C=C1/C2;
			return 1;
		}
		else if(T1&&!T2)
		{
			*T=1;
			*S=S1;
			*B=B1;
			if(C2==0)mathe=1;
			*A=A1/C2;
			*C=C1/C2;
			return 2;
		}
		else if(!T1&&T2)
		{
			if(fabs(C1)<1e-8)
			{
				*A=*T=0;
				*S="0";
				*B=1;
				*C=0;
				return 3;
			}
			else if(fabs(C2)<1e-8)
			{
				*T=1;
				*A=C1/A2;
				*B=1;
				*C=0;
				*S="/"+turnfuntostring(1,1,S2,B2,0,1);
				return 4;
			}
			else
			{
				*B=*T=1;
				*A=C1;
				*C=0;
				*S="/"+turnfuntostring(1,A2,S2,B2,C2,1);
				return 5;
			}
		}
		else
		{
			if(fabs(C1)<1e-8&&fabs(C2)<1e-8&&(S1==S2||appequal(S1,S2)))
			{
				if(fabs(B1-B2)<1e-8)
				{
					*A=*T=0;
					*S="0";
					*B=1;
					if(A2==0)mathe=0;
					*C=A1/A2;
					return 6;
				}
				else
				{
					*T=1;
					*C=0;
					if(A2==0)mathe=0;
					*A=A1/A2;
					*S=S1;
					*B=B1-B2;
					return 7;
				}
			}
			else if((S1==S2||appequal(S1,S2))&&fabs(B1-B2)<1e-8&&fabs(A1-A2)<1e-8&&fabs(C1-C2)<1e-8)
			{
				*A=*T=0;
				*C=1;
				*B=1;
				*S="0";
				return 8;
			}
			else if(fabs(C1)<1e-8&&fabs(C2)<1e-8)
			{
				*B=*T=1;
				*C=0;
				if(A2==0)mathe=0;
				*A=A1/A2;
				string tempst=turnfuntostring(1,1,S1,B1,C1,1)+"/"+turnfuntostring(1,1,S2,B2,C2,1);
				if(S1[0]=='/')
				{
					tempst=S1.substr(1,S1.length()-1);
					double tempa,tempb,tempc;
					bool tempt;
					string tempst2;
					countdiff(1,1,tempst,B1,0,1,1,S2,B2,0,&tempt,&tempa,&tempst2,&tempb,&tempc,3);
					tempst=turnfuntostring(tempt,tempa,tempst2,tempb,tempc,2);
					*S="/"+tempst;
				}
				else *S=turnfuntostring(1,1,S1,B1,C1,1)+"/"+turnfuntostring(1,1,S2,B2,C2,1);
				return 9;
			}
			else
			{
				*A=*B=*T=1;
				*C=0;
				*S=turnfuntostring(1,A1,S1,B1,C1,1)+"/"+turnfuntostring(1,A2,S2,B2,C2,1);
				return 9;
			}
		}
	}
	if(z==5)
	{
		if(!T1&&!T2)
		{
			*A=*T=0;
			*S="0";
			*B=1;
			*C=__pow(C1,C2,&mathe);
			return 1;
		}
		else if(T1&&!T2)
		{
			if(fabs(C2)<1e-8)
			{
				*A=*T=0;
				*S="0";
				*B=1;
				*C=1;
				return 2;
			}
			else if(fabs(C2-1)<1e-8)
			{
				*A=A1;
				*T=1;
				*B=B1;
				*C=C1;
				*S=S1;
				return 3;
			}
			else if(fabs(C1)<1e-8&&fabs(C2-int(C2))<1e-8&&int(C2)%2==0&&issq(S1))
			{
				int l=S1.length();
				*S=S1.substr(4,l-4);
				*A=__pow(A1,C2,&mathe);
				*T=1;
				*C=0;
				*B=B1*(int(C2)/2);
				return 10;
			}
			else if(fabs(C1)<1e-8&&fabs(C2-int(C2))<1e-8&&int(C2)%2==0&&S1[0]=='/'&&issq(S1.substr(1,S1.length()-1)))
			{
				int l=S1.length();
				*S="/"+S1.substr(5,l-5);
				*A=__pow(A1,C2,&mathe);
				*T=1;
				*C=0;
				*B=B1*(int(C2)/2);
				return 9;
			}
			else if(fabs(C1)<1e-8&&isint(C2))
			{
				*A=__pow(A1,C2,&mathe);
				*T=1;
				*B=B1*C2;
				*C=0;
				*S=S1;
				return 4;
			}
			else
			{
				*A=*T=1;
				*C=0;
				*B=C2;
				*S=turnfuntostring(1,A1,S1,B1,C1,1);
				return 5;
			}
		}
		else if(!T1&&fabs(C1)<1e-8)
		{
			*A=*T=0;
			*S="0";
			*B=1;
			*C=0;
			return 6;
		}
		else if(!T1&&fabs(C1-1)<1e-8)
		{
			*A=*T=0;
			*S="0";
			*B=*C=1;
			return 7;
		}
		else
		{
			*A=*B=*T=1;
			*C=0;
			*S=turnfuntostring(T1,A1,S1,B1,C1,3)+"^"+turnfuntostring(T2,A2,S2,B2,C2,3);
			if(!T1)return 8;
			else return 9;
		}
	}
	if(z==12||z==19)
	{
		if(!T1)
		{
			*A=*T=0;
			*B=1;
			*S="0";
			if(C1<=0)mathe=1;
			*C=log(fabs(C1));
			if(z==19)*C/=log(10);
			return 1;
		}
		else
		{
			if(fabs(C1)<1e-8)
			{
				if(fabs(A1-1)<1e-8)
				{
					*A=B1;
					if(z==19)*A/=log(10);
					*B=*T=1;
					*C=0;
					*S="ln"+turnfuntostring(1,1,S1,1,0,2);
					return 2;
				}
				else
				{
					*A=*B=*T=1;
					if(z==19)*A=1/log(10);
					*C=0;
					*S="ln"+turnfuntostring(1,A1,S1,B1,0,2);
					return 3;
				}
			}
			else
			{
				*A=*B=*T=1;
				if(z==19)*A=1/log(10);
				*C=0;
				*S="ln"+turnfuntostring(1,A1,S1,B1,C1,2);
				return 4;
			}
		}
	}
	if(z==15)
	{
		if(!T1)
		{
			*A=*T=0;
			*B=1;
			*S="0";
			if(C1<0)mathe=1;
			*C=sqrt(fabs(C1));
			return 1;
		}
		else
		{
			*A=*B=*T=1;
			*C=0;
			*S="sqrt"+turnfuntostring(1,A1,S1,B1,C1,2);
			return 2;
		}
	}
	if(z==16||z==17||z==18||(z>=21&&z<=27)||z==83||z==57||z==58||z==59)
	{
		if(!T1)
		{
			*A=*T=0;
			*B=1;
			*S="0";
			*C=count(C1,NULL,z,&mathe,0);
			return 1;
		}
		else
		{
			*A=*B=*T=1;
			*C=0;
			if(z==16)*S="sin"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==17)*S="cos"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==18)*S="tan"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==21)*S="ctg"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==22)*S="asin"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==23)*S="acos"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==24)*S="atan"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==25)*S="sh"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==26)*S="ch"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==27)*S="th"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==57)*S="ash"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==58)*S="ach"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==59)*S="ath"+turnfuntostring(1,A1,S1,B1,C1,2);
			else if(z==83)*S="actg"+turnfuntostring(1,A1,S1,B1,C1,2);
			return 2;
		}
	}
	if(z==55||z==56)
	{
		if(!T1)
		{
			*A=*T=0;
			*B=1;
			*S="0";
			*C=count(C1,NULL,z,&mathe,0);
			return 1;
		}
		else
		{
			string tempst;
			double tempa,tempb,tempc;
			bool tempt;
			countdiff(1,A1,S1,B1,C1,0,0,"0",1,0,&tempt,&tempa,&tempst,&tempb,&tempc,72-z);
			countdiff(0,0,"0",0,1,tempt,tempa,tempst,tempb,tempc,T,A,S,B,C,4);
			return 2;
		}
	}
	if(z==20)return countdiff(0,0,"0",1,Ee,T1,A1,S1,B1,C1,T,A,S,B,C,5);
	return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
void diffdebugoutput()
{
	cout<<syntaxe<<" "<<mathe<<endl;
	int i;
	for(i=1;i<=diffsum;i++)
	{
		if(diff[i].cla!=1)
		{
			if(diff[i].ope==1)cout<<" + ";
			else if(diff[i].ope==2)cout<<" - ";
			else if(diff[i].ope==3)cout<<" * ";
			else if(diff[i].ope==4)cout<<" / ";
			else if(diff[i].ope==5)cout<<" ^ ";
			else if(diff[i].ope==9)cout<<"+";
			else if(diff[i].ope==10)cout<<"-";
			else if(diff[i].ope==11)cout<<" ( ";
			else if(diff[i].ope==12)cout<<" ln( ";
			else if(diff[i].ope==14)cout<<" ) ";
			else if(diff[i].ope==15)cout<<" sqrt( ";
			else if(diff[i].ope==16)cout<<" sin( ";
			else if(diff[i].ope==17)cout<<" cos( ";
			else if(diff[i].ope==41)cout<<" ^2";
			else if(diff[i].ope==42)cout<<" ^3";
			cout<<endl;
		}
		else
		{
			if(diff[i].n[1])cout<<diff[i].fac[1]<<"*"<<diff[i].f[1]<<"^"<<diff[i].nf[1]<<" + ";
			cout<<diff[i].num[1]<<" ; ";
			if(diff[i].n[2])cout<<diff[i].fac[2]<<"*"<<diff[i].f[2]<<"^"<<diff[i].nf[2]<<" + ";
			cout<<diff[i].num[2]<<endl;
		}
	}
}
/*----------------------------------------------------------------------------------------------------------*/
void maindiff(bool *fracoutmode)
{
	char fx[1000];
	cout<<setprecision(14)<<"输入一个关于X的函数:\n";
	do
	{
		cin.getline(fx,1000,'\n');
	}while(strlen(fx)==0);
	transf(fx);
	int i,k,wtime=0;
	bool t;
	while(diffleft>1&&wtime<400&&!syntaxe&&!mathe)
	{
		wtime++;
		leftb=rightb=0;
		k=11;
		for(i=1;i<=diffsum;i++)
		{
			if(diff[i].id==1&&diff[i].cla==5)
			{
				k=diff[i].ope;
				leftb=i;
			}
			else if(diff[i].id==1&&diff[i].cla==6)
			{
				rightb=i;
				break;
			}
		}
		diff[leftb].id=diff[rightb].id=0;
		if(leftb>0&&rightb>0)diffleft-=2;
		if(leftb==0&&rightb==0)rightb=diffsum+1;
		differentiate(leftb,rightb,k);
	}
	if(syntaxe)
	{
		cout<<"表达式语法错误！\n";
		return;
	}
	if(mathe)
	{
		cout<<"数学错误！\n";
		return;
	}
	t=0;
	for(i=1;i<=diffsum;i++)
		if(diff[i].id==1&&diff[i].cla==1)
		{
			cout<<"导函数: \n";
			diffout(turnfuntostring(diff[i].n[2],diff[i].fac[2],diff[i].f[2],diff[i].nf[2],diff[i].num[2],0),fracoutmode);
			return;
		}
	if(!t)cout<<"表达式语法错误！\n";
}