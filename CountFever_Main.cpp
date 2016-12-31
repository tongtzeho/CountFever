//作者：北京大学信息科学技术学院 唐子豪
//编译环境：Microsoft Visual C++ 6.0
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <string>
#include <float.h>
#include <ctime>
#include <windows.h>

#include "CountFever_Function.h"
#include "CountFever_Input&Output.h"
#include "CountFever_SpecialMode.h"
#include "CountFever_Differentiate.h"
#include "CountFever_Plot.h"

using namespace std;

struct countfeverexpressiondata
{
	int id,cla;  //cla: 1数字 2双目 3后置 4正负号 5左括号 6右括号 7逗号 8求和变量 9X 10Y 11Z 12W 13前置 14比较 15? 16矩阵/向量/复数 17日期
	int ope;  //ope: 1+ 2- 3* 4/ 5^ 6div 7mod 8! 9+ 10- 11( 12ln( 13abs( 14) 15sq( sqrt( 16sin( 17cos( 18tan( 19lg( 20exp( 21ctg( 22asin( 23acos( 24atan( 25sh( 26ch( 27th( 28log( , 29C( 30P( 31gcd( 32lcm( 33, 34[ 35{ 36ran( 37sum(k= 38pro(k= 39I( 40% 41@ 42# 43max( 44min( 45& 46| 47xor 48_ 49< 50> 51== 52<= 53>= 54<> 55sec( 56csc( 57ash 58ach 59ath 60? 61dX) 62dY) 63dZ) 64dW) 65f( 66g( 67h( 68F( 69G( 70H( 71fib( 72sgn( 73root( 74dis( 75b& 76b| 77b^ 78~ 79<< 80>> 81rnd( 82gpa( 83actg( 84sum( 85sum2( 86on( 87on-1( 88avg( 89mid( 90n( 91avg2( 92ga( 93ha( 94pro( 95xnor 96diff( 97det( 98M( 99; 100: 101' 102" 103` 104rank( 105diag( 106tr( 107$ 108dr 109dc 110size( 111gr 112gc 113a( arg( 114re( 115im( 116sumabs( 117solve( 118floor( 119ceil( 120c( 121\ 122r( 123up( 124down( 125q 126codrp( 127sumxy( 128fc( 129A[
	int fac,r,c; //year month day
	double num,*mat;
}
data[3000];

struct countfeverfunctionvardata
{
	double start,temp,end,step;
	bool use;
	short base;
}
var[5],depth[20];

struct countfeverfunctionvaluedata
{
	double x[5],f;
}
fopti[2],bestsolve[700];

struct countfevermemorydata
{
	short fac,r,c;
	double mat[2001];
}
memans,mem[27];

int datasum,dataleft,leftbrac,rightbrac,expressionsum=0,expressionorder=0,memope,mempos,idivideorder=0,iterative=0,eqlt,solvesum=0,solvenum2;
const int idividepart[7]={2000000,160000,500,80,32,16,8};
bool funmode,ondisplay,ordinfunmode,matherror,syntaxerror,solvemode,solve1stnum,specialmode,degreemode=0,quit=0,ismatrix=0,fracoutmode=1,trianmode=0,funplot,plotting,ordinoy;
double solvelist[500]={0},randomnumlist[100]={0},*mempoint[27],codoy;
char input[5000],memexpinput[51][5000],ordinfun[7][5000],memstr[10][5000];
/*----------------------------------------------------------------------------------------------------------*/
void initializeandreadmemory(); //初始化窗口并从程序目录读取存储数据
void inputexpression(); //输入表达式
void randominitialize(); //随机初始化randomnumlist
void uppervar(); //将表达式的x,y,z,w变成大写
void clear(); //将所有存储数据清零
void showmemory(); //将所有存储数据列举出来
void dealstring(); //提取表达式中的STRx字串
bool solveon(); //判断是否启动解方程模式
void translate(); //翻译表达式，翻译结果存在结构体数组data中
void addrightbrac(); //将左括号和右括号的数目配平
void debugoutput(); //调试程序时的输出
int searchbrac(); //寻找左括号和匹配右括号
void calcexpression(const int x,const int z); //求和、求积、计算数值微积分
int maincalc(); //总引擎：每一次计算先找左右括号，再根据情况调用具体计算的函数，直到只剩一个数
void calculate(const int x,const int y,const int z); //在左右括号间计算(z为优先级)
void tranfunvar(const int x,const int y); //启动函数模式时将函数变量变成具体数值
void calcvector(const int x,const int y); //生成并计算集合/向量/数组等（最值、最大公约数、最小公倍数、求和、求平均数、求中位数、求标准差等）
void calcmatrix(const int x,const int y); //转换成矩阵类型
void calcdiag(const int x,const int y); //转换为对角矩阵
void detcm(const int x,const int n,const double number[ ]); //计算复系数方阵的行列式
void createcomplex(const int x,const double a,const double b); //生成复数
inline void matherr(const int x); //出错时的快速处理
void transpose(const int x); //矩阵转置
void adjmat(const int x); //求伴随矩阵
void invmat(const int x); //求逆矩阵
inline void oneomatrix(const int x); //将一阶矩阵变成数
void calcordinfun(const int x,const int y); //计算自定义函数
void funvarnext(); //函数模式时将函数值变成下一个
void fundatacopy(const double f,const int z); //函数模式时计算最值
bool funend(); //函数模式时判断是否结束
bool include(const int x,const int z); //判断某种运算符号是否符合优先级
void countmatrix(const int x,const int y,const int z); //对矩阵、复数运算实行计算
void complexpower(const int x,const int y); //计算复数乘方
void matrixpower(const int x,double y); //计算矩阵乘方
inline int deepen(const double x,int depnum); //解方程模式中精度向深入一层
inline int jump(int depnum); //解方程模式中精度向浅出一层
bool mulroot(const int n,const int m,const double x,const double y); //解方程模式中判断是否是重根
void newroot(const double x,const double y); //解方程模式中生成一个新的根
int solveout(int n,const double x,const double y); //解方程模式时输出方程的解
void saveanswer(); //将答案记入存储数据
void writememory(); //将存储数据写入程序目录
/*----------------------------------------------------------------------------------------------------------*/
void initializeandreadmemory()
{
	int i,j;
	SetConsoleTitle("Count Fever 1.08 中文版");
	HANDLE hOut;
	hOut=GetStdHandle(STD_OUTPUT_HANDLE);
	COORD size={1000,9999};
	SetConsoleScreenBufferSize(hOut,size);
	SMALL_RECT rc={0,0,127,31};
    SetConsoleWindowInfo(hOut,true,&rc);
	srand((unsigned)time(0));
	clear();
	mempoint[0]=&memans.mat[1];
	for(i=1;i<=26;i++)
		mempoint[i]=&mem[i].mat[1];
	ifstream infile;
	infile.open("CountFever.sav",ios::in);
	if(!infile)
	{
		cerr<<"无法读取存储数据，创建了一个新的存档\n";
		writememory();
		return;
	}
	infile>>memans.fac>>memans.r>>memans.c;
	for(i=1;i<=memans.r*memans.c&&i<=2000;i++)
		infile>>memans.mat[i];
	for(i=1;i<=26;i++)
	{
		infile>>mem[i].fac>>mem[i].r>>mem[i].c;
		for(j=1;j<=mem[i].r*mem[i].c;j++)
			infile>>mem[i].mat[j];
	}
	infile.get();
	for(i=1;i<=6;i++)
		infile.getline(ordinfun[i],5000,'\n');
	for(i=0;i<=9;i++)
		infile.getline(memstr[i],5000,'\n');
	infile.close();
	cout<<"成功读取存储数据\n";
}
/*----------------------------------------------------------------------------------------------------------*/
void inputexpression()
{
	int i,j,k,length,condisum=0,templength,quemark=0,eqlmark=0;
	char conditemp[800];
	var[1].use=var[2].use=var[3].use=var[4].use=0;
	ordinfunmode=solvemode=specialmode=trianmode=funplot=plotting=0;
	if(expressionorder==expressionsum)
	{
		expressionorder=expressionsum=0;
		for(i=0;i<5000;i++)
		    input[i]='\0';
	    cout<<endl<<"输入一条表达式或一条指令: "<<endl;
		do
		{
			cin.getline(input,5000,'\n');
			length=strlen(input);
			if(input[0]=='C'&&input[1]=='L'&&input[2]=='E'&&input[3]=='A'&&input[4]=='R')
			{
				length=0;
				clear();
				writememory();
				cout<<"成功清零存储数据"<<endl<<endl<<"输入一条表达式或一条指令: "<<endl;
			}
			else if(input[0]=='S'&&input[1]=='H'&&input[2]=='O'&&input[3]=='W'&&input[4]=='M')
			{
				length=0;
				showmemory();
				cout<<"\n输入一条表达式或一条指令: "<<endl;
			}
			else if(input[0]=='D'&&input[1]=='E'&&input[2]=='G')
			{
				length=0;
				degreemode=1;
				cout<<"使用角度制"<<endl<<endl<<"输入一条表达式或一条指令: "<<endl;
			}
			else if(input[0]=='R'&&input[1]=='A'&&input[2]=='D')
			{
				length=0;
				degreemode=0;
				cout<<"使用弧度制"<<endl<<endl<<"输入一条表达式或一条指令: "<<endl;
			}
			else if(input[0]=='D'&&input[1]=='E'&&input[2]=='C')
			{
				length=0;
				fracoutmode=0;
				cout<<"输出小数"<<endl<<endl<<"输入一条表达式或一条指令: "<<endl;
			}
			else if(input[0]=='F'&&input[1]=='R'&&input[2]=='A')
			{
				length=0;
				fracoutmode=1;
				cout<<"输出分数或有理根式"<<endl<<endl<<"输入一条表达式或一条指令: "<<endl;
			}
			else if(input[0]=='C'&&input[1]=='L'&&input[2]=='S'&&input[3]=='C'&&input[4]=='R')
			{
				length=0;
				system("cls");
				cout<<"成功清屏\n\n输入一条表达式或一条指令: "<<endl;
			}
		}while(length==0);
		if(input[0]=='E'&&input[1]=='X'&&input[2]=='P'&&input[3]=='R'&&input[4]=='E')
		{
			cout<<"表达式的数量? ( 1 -- 50 )       （缺省值为5）\n";
			expressionsum=int(inputvardata(mempoint,5,1));
			if(expressionsum>50||expressionsum<=0)expressionsum=0;
			for(i=1;i<=expressionsum;i++)
			{
				cout<<"输入第"<<i<<"条表达式 :\n";
				for(j=0;j<5000;j++)
					memexpinput[i][j]='\0';
				do
				{
					cin.getline(memexpinput[i],5000,'\n');
				}while(strlen(memexpinput[i])==0);
			}
		}
	}
	if(expressionsum>0)
	{
		expressionorder++;
		cout<<endl<<"表达式"<<expressionorder<<" :\n";
		for(i=0;i<5000;i++)
			input[i]=memexpinput[expressionorder][i];
		length=strlen(input);
		for(i=0;i<length;i++)
			cout<<input[i];
		cout<<endl;
	}
	dealstring();
	uppervar();
	length=strlen(input);
	if((input[0]=='S'&&input[1]=='O'&&input[2]=='L'&&input[3]=='V'&&input[4]=='E')||solveon())
	{
		solvemode=funmode=1;
		if(input[0]=='S')
		{
			cout<<"输入一条关于X的方程: \n";
		    for(i=0;i<5000;i++)
			    input[i]='\0';
		    do
			{
				cin.getline(input,5000,'\n');
				dealstring();
			}while(strlen(input)==0&&!solveon());
		}
		else dealstring();
		length=strlen(input);
		for(i=0;i<length;i++)
			if(input[i]=='='&&input[i+1]!='='&&input[i-1]!='='&&input[i-1]!='k'&&input[i-1]!='>'&&input[i-1]!='<')break;
		for(j=length-1;j>=0;j--)
		{
			if(j>i)input[j+3]=input[j];
			else if(j==i)input[j+2]='-';
			else if(j<i)input[j+1]=input[j];
		}
		input[0]=input[i+3]='q';
		input[length+3]=input[i+1]=')';
		length=strlen(input);
		var[0].use=var[1].use=1;
		cout<<"X≥?       （缺省值为-20）\n";
		var[1].start=var[1].temp=inputvardata(mempoint,-20,1);
		cout<<"X≤?       （缺省值为20）\n";
		var[1].end=inputvardata(mempoint,20,1);
		while(var[1].end<=var[1].start)
		{
			cout<<"出错！重新输入 :       （缺省值为20）\n";
			var[1].end=inputvardata(mempoint,20,1);
		}
		solvesum=99;
		var[1].step=(var[1].end-var[1].start)/5000;
		var[1].end+=(var[1].step/100);
		depth[0]=var[1];
		for(i=1;i<10;i++)
			depth[i].start=depth[i].end=depth[i].step=depth[i].temp=0;
		uppervar();
		return;
	}
	else if(input[0]=='P'&&input[1]=='R'&&input[2]=='I'&&input[3]=='F'&&input[4]=='A')
	{
		specialmode=1;
		primefactorize();
		return;
	}
	else if(input[0]=='T'&&input[1]=='R'&&input[2]=='I'&&input[3]=='A'&&input[4]=='N')
	{
		specialmode=trianmode=1;
		triangle(degreemode,fracoutmode,mempoint);
		return;
	}
	else if(input[0]=='H'&&input[1]=='E'&&input[2]=='X'&&input[3]=='C'&&input[4]=='O')
	{
		specialmode=1;
		hexconver();
		return;
	}
	else if(input[0]=='D'&&input[1]=='I'&&input[2]=='F'&&input[3]=='F')
	{
		specialmode=1;
		maindiff(&fracoutmode);
		return;
	}
	else if(input[0]=='H'&&input[1]=='E'&&input[2]=='L'&&input[3]=='P')
	{
		specialmode=1;
		help();
	}
	else if((input[0]=='E'&&input[1]=='X'&&input[2]=='I'&&input[3]=='T')||(input[0]=='O'&&input[1]=='F'&&input[2]=='F')||(input[0]=='Q'&&input[1]=='U'&&input[2]=='I'&&input[3]=='T'))
	{
		quit=1;
		return;
	}
	if(((input[0]>='f'&&input[0]<='h')||(input[0]>='F'&&input[0]<='H'))&&input[1]=='('&&input[2]>=87&&input[2]<=90&&(input[3]==','||input[3]==')'))
	{
		funmode=1;
		if(input[4]=='=')i=5;
	    else if(input[6]=='=')i=7;
	    else if(input[8]=='=')i=9;
	    else if(input[10]=='=')i=11;
		else i=10;
		eqlmark=i==10?0:i;
		if(eqlmark>0)
		{
		    k=input[0]<91?input[0]-66:input[0]-101;
		    for(j=0;j<5000;j++)
			    ordinfun[k][j]='\0';
		    for(j=0;j<length;j++)
			    ordinfun[k][j]=input[j];
			writememory();
		}
		if(input[length-1]==92)return;
	    for(j=2;j<i;j++)
		{
		    if(input[j]=='X')var[1].use=1;
		    else if(input[j]=='Y')var[2].use=1;
		    else if(input[j]=='Z')var[3].use=1;
		    else if(input[j]=='W')var[4].use=1;
		}
		if(eqlmark==0)
		{
			i=j=0;
			while(input[i]!=')')
			{
				if(!((input[i]>=70&&input[i]<=72)||(input[i]>=102&&input[i]<=104)||(input[i]>='W'&&input[i]<='Z')||input[i]==','||input[i]=='('))j=1;
				i++;
				if(i==5000)break;
			}
			i=0;
			if(j==0)ordinfunmode=1;
			else return;
		}
		for(j=1;j<=4;j++)
			if(var[j].use==1)
			{
				cout<<char(j<=3?j+87:87)<<"的起始值?      （缺省值为-10）\n";
				var[j].start=inputvardata(mempoint,-10,1);
				cout<<char(j<=3?j+87:87)<<"的终止值?      （缺省值为10）\n";
				var[j].end=inputvardata(mempoint,10,1);
				cout<<char(j<=3?j+87:87)<<"的间隔?        （缺省值为1）\n";
				var[j].step=inputvardata(mempoint,1,1);
				while(var[j].step<=0)
				{
					cout<<"出错！间隔必须大于零！重新输入 :     （缺省值为1）\n";
					var[j].step=inputvardata(mempoint,1,1);
				}
				var[j].temp=var[j].start;
				var[j].end+=1e-9;
			}
		cout<<"限定条件的数量? ( 0 -- 20 )    （缺省值为0）\n";
		cin.getline(conditemp,800,'\n');
		condisum=atoi(conditemp);
		if(strlen(conditemp)==0||abs(condisum)>20)condisum=0;
		for(j=1;j<=condisum;j++)
		{
			cout<<"输入第"<<j<<"个限定条件 :\n";
			cin.getline(conditemp,800,'\n');
			templength=strlen(conditemp);
			for(k=length-1;k>=i;k--)
				input[k+templength+3]=input[k];
			if(j>1)input[i-1]='&';
			input[i]='(';
			input[i+templength+1]=')';
			input[i+templength+2]='?';
			for(k=1;k<=templength;k++)
				input[i+k]=conditemp[k-1];
			i+=(templength+3);
			length+=(templength+3);
		}
		if(condisum==0&&var[1].use+var[2].use+var[3].use+var[4].use==1)
		{
			cout<<"是否绘制函数图象? (y/n)     （默认为是）\n";
			cin.getline(conditemp,800,'\n');
			if(!(conditemp[0]=='N'||conditemp[0]=='n'))funplot=1;
			if(funplot)
			{
				cout<<"坐标原点的纵坐标?      （缺省值为系统设定值）\n";
				codoy=inputvardata(mempoint,1.2345e+308,1);
				ordinoy=codoy!=1.2345e+308;
			}
		}
	}
	for(i=0;i<length;i++)
		if(input[i]=='?')
		{
			quemark=i;
			break;
		}
	if(quemark)
	{
		for(i=length-1;i>=quemark;i--)
			input[i+2]=input[i];
		input[quemark+1]=')';
		for(i=quemark-1;i>=eqlmark;i--)
			input[i+1]=input[i];
		input[eqlmark]='(';
	}
	length=strlen(input);
	uppervar();
}
/*----------------------------------------------------------------------------------------------------------*/
void randominitialize()
{
	int i;
	for(i=0;i<100;i++)
		randomnumlist[i]=ran(9580032000.0);
}
/*----------------------------------------------------------------------------------------------------------*/
void uppervar()
{
	int i,length=strlen(input);
	for(i=0;i<length;i++)
	{
		if(input[i]=='y'||input[i]=='z'||input[i]=='w')input[i]-=32;
		else if(input[i]=='x'&&input[i-1]!='a'&&input[i+1]!='p'&&input[i+1]!='o'&&input[i+1]!='n')input[i]='X';
	}
}
/*----------------------------------------------------------------------------------------------------------*/
void clear()
{
	const char ord[ ]=" (X)=0\0";
	int i,j,length=strlen(ord);
	memans.fac=0;
	memans.r=memans.c=1;
	for(i=0;i<=2000;i++)
		memans.mat[i]=0;
	for(i=0;i<=26;i++)
	{
		mem[i].fac=0;
		mem[i].r=mem[i].c=1;
		for(j=0;j<=2000;j++)
			mem[i].mat[j]=0;
	}
	mempos=memope=0;
	for(i=1;i<=6;i++)
		for(j=0;j<5000;j++)
			ordinfun[i][j]='\0';
	ordinfun[1][0]='f';
	ordinfun[2][0]='g';
	ordinfun[3][0]='h';
	ordinfun[4][0]='F';
	ordinfun[5][0]='G';
	ordinfun[6][0]='H';
	for(i=1;i<=6;i++)
		for(j=1;j<length;j++)
			ordinfun[i][j]=ord[j];
	for(i=0;i<=9;i++)
		for(j=0;j<5000;j++)
			memstr[i][j]='\0';
	for(i=0;i<=9;i++)
		memstr[i][0]='0';
}
/*----------------------------------------------------------------------------------------------------------*/
void showmemory()
{
	ifstream infile;
	infile.open("CountFever.sav",ios::in);
	if(!infile)
	{
		cerr<<"出错！无法读取存储数据！\n";
		return;
	}
	double temp,temp2;
	int r,c,fac,i,j,k;
	for(i=0;i<=26;i++)
	{
	    infile>>fac>>r>>c;
		if(i==0)cout<<"ANS=";
		else cout<<"M"<<char(i+64)<<"=";
	    if(fac==0)
		{
		    infile>>temp;
		    cout<<setprecision(12)<<temp<<endl;
		}
	    else if(fac==1&&r==1&&c==2)
		{
		    infile>>temp>>temp2;
		    complexdecout(temp,temp2);
		    cout<<endl;
		}
	    else
		{
		    if(c>70)
			{
				cout<<"(矩阵太大无法输出)"<<endl;
				for(j=1;j<=r*c;j++)
					infile>>temp;
			}
			else
			{
				cout<<endl;
		        for(j=1;j<=r;j++)
				{
			        if(r==1)cout<<"（";
			        else if(j==1)cout<<"┌";
			        else if(j==r)cout<<"└";
			        else cout<<"│";
			        for(k=1;k<=c;k++)
					{
				        infile>>temp;
				        cout<<setprecision(6)<<setw(14)<<temp;
					}
			        if(r==1)cout<<"    ）";
			        else if(j==1)cout<<"    ┐";
			        else if(j==r)cout<<"    ┘";
			        else cout<<"    │";
			        if(j==r)cout<<r<<"×"<<c<<endl;
			        else cout<<endl;
				}
			}
		}
	}
	infile.get();
	char str[5000];
	i=0;
	while(infile.getline(str,5000,'\n'))
	{
		i++;
		if(i<=6)cout<<str<<endl;
		else cout<<"STR"<<i-7<<":"<<str<<endl;
	}
	cout<<"成功列举存储数据\n";
	infile.close();
}
/*----------------------------------------------------------------------------------------------------------*/
void dealstring()
{
	bool t=0;
	int i,j,length;
	length=strlen(input);
	for(i=0;i<length;i++)
		if(input[i]=='S'&&input[i+1]=='T'&&input[i+2]=='R'&&input[i+3]>='0'&&input[i+3]<='9'&&(i<=1||(i>=2&&input[i-1]!='>'&&input[i-2]!='-')))
		{
			t=1;
			break;
		}
	ofstream outfile;
	ifstream infile;
	if(t)
	{
		outfile.open("temp.txt",ios::out);
		for(i=0;i<length; )
		{
		    if(input[i]=='S'&&input[i+1]=='T'&&input[i+2]=='R'&&input[i+3]>='0'&&input[i+3]<='9'&&(i<=1||(i>=2&&input[i-1]!='>'&&input[i-2]!='-')))
			{
				for(j=0;j<strlen(memstr[input[i+3]-48]);j++)
					outfile<<memstr[input[i+3]-48][j];
				i+=4;
			}
			else
			{
				outfile<<input[i];
				i++;
			}
		}
		for(i=0;i<5000;i++)
			input[i]='\0';
		outfile.close();
		infile.open("temp.txt",ios::in);
		infile.getline(input,5000,'\n');
		infile.close();
		length=strlen(input);
		for(i=0;i<length;i++)
			cout<<input[i];
		cout<<endl;
		remove("temp.txt");
	}
	if(length>=6&&input[length-1]>=48&&input[length-1]<=57&&input[length-2]=='R'&&input[length-3]=='T'&&input[length-4]=='S'&&input[length-5]=='>'&&input[length-6]=='-')
	{
		j=input[length-1]-48;
		for(i=length-6;i<=length;i++)
			input[i]='\0';
		length=strlen(input);
		for(i=0;i<5000;i++)
			memstr[j][i]='\0';
		for(i=0;i<length;i++)
			memstr[j][i]=input[i];
	}
}
/*----------------------------------------------------------------------------------------------------------*/
bool solveon()
{
	int i,length=strlen(input),eqlsum=0;
	if(((input[0]>='f'&&input[0]<='h')||(input[0]>='F'&&input[0]<='H'))&&input[1]=='('&&input[2]>=87&&input[2]<=90&&(input[3]==','||input[3]==')'))return false;
	for(i=0;i<length;i++)
		if(input[i]=='='&&input[i+1]!='='&&input[i-1]!='='&&input[i-1]!='k'&&input[i-1]!='>'&&input[i-1]!='<')eqlsum++;
	if(eqlsum==1)return true;
	return false;
}
/*----------------------------------------------------------------------------------------------------------*/
void translate()
{
	const double ER=0.57721566490153286061; //R
    const double Kk=2.6854520010653064453; //K
    const double CWY=1.3035772690342963913; //CWY
    const double Cc=299792458.0; //CO
    const double Gg=6.6726e-11; //G
    const double ee=1.6021892e-19; //E
    const double NA=6.02213674e+23; //NA
    const double Hh=6.626069e-34; //H
    const double me=9.109939e-31; //me
    const double mp=1.6726485e-27; //mp
    const double mE=5.997e+24; //mE
    const double RD=1.096776e+7; //RY
    const double Ro=8.31441; //RO
    const double ATM=101325.0; //ATM
    const double gg=9.80665; //g
	const double E0=8.854187818e-12; //e0
	datasum=dataleft=mempos=funmode=0;
	int i=0,j,length,t=0,t2=0,templength,pointsum;
	char temp[500];
	length=strlen(input);
	if(solvemode==1)funmode=1;
	if(syntaxerror)return;
	else if(((input[0]>='f'&&input[0]<='h')||(input[0]>='F'&&input[0]<='H'))&&input[1]=='('&&input[2]>=87&&input[2]<=90&&(input[3]==','||input[3]==')'))
	{
		funmode=1;
		if(input[4]=='=')i=5;
		else if(input[6]=='=')i=7;
		else if(input[8]=='=')i=9;
		else if(input[10]=='=')i=11;
		else i=funmode=0;
	}
	if(ordinfunmode==1)funmode=1;
	if(input[length-1]==92)
	{
		length--;
		input[length]='\0';
		funmode=0;
	}
	while(i<length)
	{
		while(input[i]==' ')
			i++;
		if(input[i]=='+'||(input[i]=='-'&&input[i+1]!='>'))
		{
			if(t==0)
			{
				datasum++;
				data[datasum].id=datasum;
				data[datasum].cla=4;
				if(input[i]=='+')data[datasum].ope=9;
				else data[datasum].ope=10;
			}
			if(t==1)
			{
				datasum++;
				data[datasum].id=datasum;
				data[datasum].cla=2;
				if(input[i]=='+')data[datasum].ope=1;
				else data[datasum].ope=2;
				t=0;
			}
			t2=0;
		}
		else if(input[i]=='*'||input[i]=='/'||input[i]=='^'||(input[i]=='d'&&input[i+1]=='i'&&input[i+2]=='v')||(input[i]=='m'&&input[i+1]=='o')||input[i]=='&'||input[i]=='|'||(input[i]=='x'&&(input[i+1]=='o'||input[i+1]=='n'))||(input[i]=='<'&&input[i+1]=='<')||(input[i]=='>'&&input[i+1]=='>')||input[i]=='b'||input[i]=='`'||input[i]=='$'||(input[i]=='d'&&(input[i+1]=='r'||input[i+1]=='c'))||(input[i]=='g'&&(input[i+1]=='r'||(input[i+1]=='c'&&input[i+2]!='d')))||input[i]==92||(input[i]=='C'&&(input[i+1]!='W'||input[i+2]!='Y')&&input[i+1]!='O')||(input[i]=='A'&&input[i+1]!='N'&&input[i+1]!='T'&&input[i+1]!='[')||input[i]=='P')
		{
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=2;
			if(solvemode&&(input[i]=='&'||input[i]=='|'||input[i]=='x'||input[i]=='b'||input[i]=='>'||input[i]=='<'||(input[i]=='='&&input[i+1]=='=')))
			{
				syntaxerror=1;
				return;
			}
			if(input[i]=='*')data[datasum].ope=3;
			else if(input[i]=='/')data[datasum].ope=4;
			else if(input[i]=='^')data[datasum].ope=5;
			else if(input[i]=='d'&&input[i+1]=='i')
			{
				data[datasum].ope=6;
				i+=2;
			}
			else if(input[i]=='m')
			{
				data[datasum].ope=7;
				i+=2;
			}
			else if(input[i]=='&')data[datasum].ope=45;
			else if(input[i]=='|')data[datasum].ope=46;
			else if(input[i]=='x'&&input[i+1]=='o')
			{
				data[datasum].ope=47;
				i+=2;
			}
			else if(input[i]=='b')
			{
				i++;
				if(input[i]=='&')data[datasum].ope=75;
				else if(input[i]=='|')data[datasum].ope=76;
				else if(input[i]=='^')data[datasum].ope=77;
			}
			else if(input[i]=='<')
			{
				data[datasum].ope=79;
				i++;
			}
			else if(input[i]=='>')
			{
				data[datasum].ope=80;
				i++;
			}
			else if(input[i]=='x'&&input[i+1]=='n')
			{
				data[datasum].ope=95;
				i+=3;
			}
			else if(input[i]=='`')data[datasum].ope=103;
			else if(input[i]=='$')data[datasum].ope=107;
			else if(input[i]=='d'&&input[i+1]=='r')
			{
				data[datasum].ope=108;
				i++;
			}
			else if(input[i]=='d'&&input[i+1]=='c')
			{
				data[datasum].ope=109;
				i++;
			}
			else if(input[i]=='g'&&input[i+1]=='r')
			{
				data[datasum].ope=111;
				i++;
			}
			else if(input[i]=='g'&&input[i+1]=='c')
			{
				data[datasum].ope=112;
				i++;
			}
			else if(input[i]==92)data[datasum].ope=121;
			else if(input[i]=='C')data[datasum].ope=29;
			else if(input[i]=='A'||input[i]=='P')data[datasum].ope=30;
			t=t2=0;
		}
		else if(input[i]=='!'||input[i]=='%'||input[i]=='@'||input[i]=='#'||input[i]==39||input[i]==34)
		{
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=3;
			if(input[i]=='!')data[datasum].ope=8;
			else if(input[i]=='%')data[datasum].ope=40;
			else if(input[i]=='@')data[datasum].ope=41;
			else if(input[i]=='#')data[datasum].ope=42;
			else if(input[i]==39)data[datasum].ope=101;
			else if(input[i]==34)data[datasum].ope=102;
			t=1;
			if(input[i]=='@'||input[i]=='#')t2=1;
			else t2=0;
		}
		else if(input[i]=='_'||input[i]=='~')
		{
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=13;
			if(input[i]=='_')data[datasum].ope=48;
			else if(input[i]=='~')data[datasum].ope=78;
			t=t2=0;
		}
		else if((input[i]=='<'&&input[i+1]!='<')||(input[i]=='>'&&input[i+1]!='>')||(input[i]=='='&&input[i+1]=='='))
		{
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=14;
			if(solvemode)
			{
				syntaxerror=1;
				return;
			}
			if(input[i]=='<'&&input[i+1]=='=')
			{
				data[datasum].ope=52;
				i++;
			}
			else if(input[i]=='>'&&input[i+1]=='=')
			{
				data[datasum].ope=53;
				i++;
			}
			else if(input[i]=='=')
			{
				data[datasum].ope=51;
				i++;
			}
			else if(input[i]=='<'&&input[i+1]=='>')
			{
				data[datasum].ope=54;
				i++;
			}
			else if(input[i]=='<')data[datasum].ope=49;
			else if(input[i]=='>')data[datasum].ope=50;
			t=t2=0;
		}
		else if(input[i]=='?')
		{
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=15;
			data[datasum].ope=60;
			t=t2=0;
		}
		else if(input[i]=='a'||input[i]=='l'||input[i]=='('||input[i]=='s'||input[i]=='c'||input[i]=='t'||(input[i]=='e'&&input[i+1]=='x')||(input[i]=='g'&&(input[i+1]=='c'||input[i+1]=='p'||input[i+1]=='a'||input[i+1]=='('))||input[i]=='['||input[i]=='{'||(input[i]=='p'&&input[i+1]=='r')||input[i]=='I'||(input[i]=='m'&&(input[i+1]=='a'||input[i+1]=='i'))||input[i]=='f'||(input[i]=='h'&&(input[i+1]=='('||input[i+1]=='a'))||input[i]=='F'||(input[i]=='G'&&input[i+1]=='(')||(input[i]=='H'&&input[i+1]=='(')||(input[i]=='r'&&(input[i+1]=='o'||input[i+1]=='e'||input[i+1]=='a'||input[i+1]=='n'||input[i+1]=='('||input[i+1]=='['))||(input[i]=='d'&&(input[i+1]=='i'||input[i+1]=='e'||input[i+1]=='o'))||(input[i]=='o'&&input[i+1]=='n')||((input[i]=='M'||input[i]=='n')&&(input[i+1]=='('||input[i+1]=='['))||(input[i]=='i'&&input[i+1]=='m')||(input[i]=='u'&&input[i+1]=='p')||(solvemode&&input[i]=='q')||(input[i]=='A'&&input[i+1]=='['))
		{
			if(t2==1)
			{
				datasum++;
				data[datasum].id=datasum;
			    data[datasum].cla=2;
				data[datasum].ope=3;
			}
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=5;
			if(input[i]=='a'&&input[i+1]=='b')
			{
				data[datasum].ope=13;
				i+=3;
			}
			else if(input[i]=='l'&&input[i+1]=='n')
			{
				data[datasum].ope=12;
				i+=2;
			}
			else if(input[i]=='s'&&input[i+1]=='q')
			{
				data[datasum].ope=15;
				if(input[i+2]=='r')i+=4;
				else if(input[i+2]=='(')i+=2;
				else
				{
					syntaxerror=1;
					return;
				}
			}
			else if(input[i]=='s'&&input[i+1]=='i'&&input[i+2]=='n')
			{
				data[datasum].ope=16;
				i+=3;
			}
			else if(input[i]=='c'&&input[i+1]=='o'&&input[i+2]=='s')
			{
				data[datasum].ope=17;
				i+=3;
			}
			else if(input[i]=='t'&&(input[i+1]=='a'||input[i+1]=='g'))
			{
				data[datasum].ope=18;
				input[i+1]=='a'?i+=3:i+=2;
			}
			else if(input[i]=='l'&&input[i+1]=='g')
			{
				data[datasum].ope=19;
				i+=2;
			}
			else if(input[i]=='e'&&input[i+1]=='x')
			{
				data[datasum].ope=20;
				i+=3;
			}
			else if(input[i]=='c'&&(input[i+1]=='t'||(input[i+1]=='o'&&input[i+2]=='t')))
			{
				data[datasum].ope=21;
				i+=3;
			}
			else if(input[i]=='a'&&input[i+1]=='s'&&input[i+2]=='i')
			{
				data[datasum].ope=22;
				i+=4;
			}
			else if(input[i]=='a'&&input[i+1]=='c'&&input[i+2]=='o'&&input[i+3]=='s')
			{
				data[datasum].ope=23;
				i+=4;
			}
			else if(input[i]=='a'&&input[i+1]=='t'&&(input[i+2]=='a'||input[i+2]=='g'))
			{
				data[datasum].ope=24;
				input[i+2]=='a'?i+=4:i+=3;
			}
			else if(input[i]=='s'&&input[i+1]=='h')
			{
				data[datasum].ope=25;
				i+=2;
			}
			else if(input[i]=='c'&&input[i+1]=='h')
			{
				data[datasum].ope=26;
				i+=2;
			}
			else if(input[i]=='t'&&input[i+1]=='h')
			{
				data[datasum].ope=27;
				i+=2;
			}
			else if(input[i]=='l'&&input[i+1]=='o')
			{
				data[datasum].ope=28;
				i+=3;
			}
			else if(input[i]=='g'&&input[i+1]=='c')
			{
				data[datasum].ope=31;
				i+=3;
			}
			else if(input[i]=='l'&&input[i+1]=='c')
			{
				data[datasum].ope=32;
				i+=3;
			}
			else if(input[i]=='[')data[datasum].ope=34;
			else if(input[i]=='{')data[datasum].ope=35;
			else if(input[i]=='r'&&input[i+1]=='a'&&(input[i+3]=='('||input[i+3]=='['))
			{
				data[datasum].ope=36;
				i+=3;
			}
			else if(input[i]=='s'&&input[i+1]=='u'&&input[i+4]=='k')
			{
				data[datasum].ope=37;
				i+=5;
			}
			else if(input[i]=='p'&&input[i+1]=='r'&&input[i+4]=='k')
			{
				data[datasum].ope=38;
				i+=5;
			}
			else if(input[i]=='I'&&(input[i+1]=='('||input[i+2]=='('))
			{
				data[datasum].ope=39;
				if(input[i+1]=='(')
				{
					data[datasum].fac=-1;
				    i++;
				}
				else if(48<=input[i+1]&&input[i+1]<=57)
				{
					data[datasum].fac=input[i+1]-48;
					if(data[datasum].fac>=7)
					{
						syntaxerror=1;
						return;
					}
					i+=2;
				}
				else
				{
					syntaxerror=1;
					return;
				}
			}
			else if(input[i]=='m'&&input[i+1]=='a')
			{
				data[datasum].ope=43;
				i+=3;
			}
			else if(input[i]=='m'&&input[i+1]=='i'&&input[i+2]=='n')
			{
				data[datasum].ope=44;
				i+=3;
			}
			else if(input[i]=='s'&&input[i+1]=='e')
			{
				data[datasum].ope=55;
				i+=3;
			}
			else if(input[i]=='c'&&input[i+1]=='s')
			{
				data[datasum].ope=56;
				i+=3;
			}
			else if(input[i]=='a'&&input[i+1]=='s')
			{
				data[datasum].ope=57;
				i+=3;
			}
			else if(input[i]=='a'&&input[i+1]=='c'&&input[i+2]=='h')
			{
				data[datasum].ope=58;
				i+=3;
			}
			else if(input[i]=='a'&&input[i+1]=='t')
			{
				data[datasum].ope=59;
				i+=3;
			}
			else if(input[i]=='f'&&input[i+1]=='(')
			{
				data[datasum].ope=65;
				i++;
			}
			else if(input[i]=='g'&&input[i+1]=='(')
			{
				data[datasum].ope=66;
				i++;
			}
			else if(input[i]=='h'&&input[i+1]=='(')
			{
				data[datasum].ope=67;
				i++;
			}
			else if(input[i]=='F'&&input[i+1]=='(')
			{
				data[datasum].ope=68;
				i++;
			}
			else if(input[i]=='G'&&input[i+1]=='(')
			{
				data[datasum].ope=69;
				i++;
			}
			else if(input[i]=='H'&&input[i+1]=='(')
			{
				data[datasum].ope=70;
				i++;
			}
			else if(input[i]=='f'&&input[i+1]=='i')
			{
				data[datasum].ope=71;
				i+=3;
			}
			else if(input[i]=='s'&&input[i+1]=='g')
			{
				data[datasum].ope=72;
				i+=3;
			}
			else if(input[i]=='r'&&input[i+1]=='o')
			{
				data[datasum].ope=73;
				i+=4;
			}
			else if(input[i]=='d'&&input[i+1]=='i'&&input[i+2]=='s')
			{
				data[datasum].ope=74;
				i+=3;
			}
			else if(input[i]=='r'&&input[i+1]=='n')
			{
				data[datasum].ope=81;
				i+=3;
			}
			else if(input[i]=='g'&&input[i+1]=='p'&&input[i+2]=='a')
			{
				data[datasum].ope=82;
				if(input[i+3]=='('||input[i+3]=='[')
				{
					i+=3;
					data[datasum].fac=0;
				}
				else
				{
					if(input[i+3]=='a')data[datasum].fac=1;
					else if(input[i+3]=='p')data[datasum].fac=2;
					else if(input[i+3]=='s')data[datasum].fac=3;
					else if(input[i+3]=='r')data[datasum].fac=4;
					else syntaxerror=1;
					i+=4;
				}
			}
			else if(input[i]=='a'&&input[i+1]=='c'&&(input[i+2]=='t'||input[i+2]=='o'))
			{
				data[datasum].ope=83;
				i+=4;
			}
			else if(input[i]=='s'&&input[i+1]=='u'&&(input[i+3]=='('||input[i+3]=='['||input[i+3]=='f'))
			{
				data[datasum].ope=84;
				if(input[i+3]=='f')
				{
					i+=4;
					data[datasum].fac=1;
				}
				else 
				{
					i+=3;
					data[datasum].fac=0;
				}
			}
			else if(input[i]=='s'&&input[i+1]=='u'&&input[i+3]>='1'&&input[i+3]<='9')
			{
				data[datasum].ope=85;
				data[datasum].fac=input[i+3]-'0';
				if(input[i+4]=='f')
				{
					i+=5;
					data[datasum].fac+=10;
				}
				else i+=4;
			}
			else if(input[i]=='o'&&input[i+1]=='n'&&(input[i+2]=='('||input[i+2]=='['||input[i+2]=='f'||input[i+2]=='X'||input[i+2]=='Y'))
			{
				data[datasum].ope=86;
				if(input[i+2]=='f'||input[i+2]=='X'||input[i+2]=='Y')
				{
					if(input[i+2]=='f')data[datasum].fac=1;
					else if(input[i+2]=='X')data[datasum].fac=2;
					else data[datasum].fac=3;
					i+=3;
				}
				else
				{
					i+=2;
					data[datasum].fac=0;
				}
			}
			else if(input[i]=='o'&&input[i+1]=='n'&&input[i+2]=='-'&&input[i+3]=='1')
			{
				data[datasum].ope=87;
				if(input[i+4]=='f'||input[i+4]=='X'||input[i+4]=='Y')
				{
					if(input[i+4]=='f')data[datasum].fac=1;
					else if(input[i+4]=='X')data[datasum].fac=2;
					else data[datasum].fac=3;
	                i+=5;
				}
				else if(input[i+4]=='('||input[i+4]=='[')
				{
					i+=4;
					data[datasum].fac=0;
				}
				else
				{
					syntaxerror=1;
					return;
				}
			}
			else if(input[i]=='a'&&input[i+1]=='v'&&(input[i+3]=='('||input[i+3]=='['||input[i+3]=='f'||input[i+3]=='X'||input[i+3]=='Y'))
			{
				data[datasum].ope=88;
				if(input[i+3]=='f'||input[i+3]=='X'||input[i+3]=='Y')
				{
					if(input[i+3]=='f')data[datasum].fac=1;
					else if(input[i+3]=='X')data[datasum].fac=2;
					else data[datasum].fac=3;
					i+=4;
				}
				else
				{
					i+=3;
					data[datasum].fac=0;
				}
			}
			else if(input[i]=='m'&&input[i+1]=='i'&&input[i+2]=='d')
			{
				data[datasum].ope=89;
				if(input[i+3]=='f')
				{
					i+=4;
					data[datasum].fac=1;
				}
				else
				{
					i+=3;
					data[datasum].fac=0;
				}
			}
			else if(input[i]=='n'&&(input[i+1]=='('||input[i+1]=='['||input[i+1]=='f'))
			{
				data[datasum].ope=90;
				if(input[i+1]=='f')
				{
					i+=2;
					data[datasum].fac=1;
				}
				else
				{
					i++;
					data[datasum].fac=0;
				}
			}
			else if(input[i]=='a'&&input[i+1]=='v'&&input[i+3]=='2')
			{
				data[datasum].ope=91;
				if(input[i+4]=='f')
				{
				    i+=5;
					data[datasum].fac=1;
				}
				else
				{
					i+=4;
					data[datasum].fac=0;
				}
			}
			else if(input[i]=='g'&&input[i+1]=='a')
			{
				data[datasum].ope=92;
				if(input[i+2]=='f')
				{
					i+=3;
					data[datasum].fac=1;
				}
				else
				{
					i+=2;
					data[datasum].fac=0;
				}
			}
			else if(input[i]=='h'&&input[i+1]=='a')
			{
				data[datasum].ope=93;
				if(input[i+2]=='f')
				{
					i+=3;
					data[datasum].fac=1;
				}
				else
				{
					i+=2;
					data[datasum].fac=0;
				}
			}
			else if(input[i]=='p'&&input[i+1]=='r'&&(input[i+3]=='('||input[i+3]=='['||input[i+3]=='f'))
			{
				data[datasum].ope=94;
				if(input[i+3]=='f')
				{
					i+=4;
					data[datasum].fac=1;
				}
				else 
				{
					i+=3;
					data[datasum].fac=0;
				}
			}
			else if(input[i]=='d'&&input[i+1]=='i'&&input[i+2]=='f')
			{
				data[datasum].ope=96;
				i+=4;
			}
			else if(input[i]=='d'&&input[i+1]=='e')
			{
				data[datasum].ope=97;
				i+=3;
			}
			else if(input[i]=='M'&&(input[i+1]=='('||input[i+1]=='['))
			{
				data[datasum].ope=98;
				i++;
			}
			else if(input[i]=='r'&&input[i+1]=='a'&&input[i+3]=='k')
			{
				data[datasum].ope=104;
				i+=4;
			}
			else if(input[i]=='d'&&input[i+1]=='i'&&input[i+2]=='a')
			{
				data[datasum].ope=105;
				i+=4;
			}
			else if(input[i]=='t'&&input[i+1]=='r')
			{
				data[datasum].ope=106;
				i+=2;
			}
			else if(input[i]=='s'&&input[i+1]=='i'&&(input[i+2]=='z'||input[i+2]=='Z'))
			{
				data[datasum].ope=110;
				i+=4;
			}
			else if(input[i]=='a'&&(input[i+1]=='('||input[i+1]=='['||input[i+1]=='r'))
			{
				data[datasum].ope=113;
				if(input[i+1]=='r')i+=3;
				else i++;
			}
			else if(input[i]=='r'&&input[i+1]=='e')
			{
				data[datasum].ope=114;
				i+=2;
			}
			else if(input[i]=='i'&&input[i+1]=='m')
			{
				data[datasum].ope=115;
				i+=2;
			}
			else if(input[i]=='s'&&input[i+1]=='u'&&input[i+4]=='b'&&input[i+5]=='s'&&(input[i+6]=='f'||input[i+6]=='['||input[i+6]=='('))
			{
				data[datasum].ope=116;
				if(input[i+6]=='f')
				{
					i+=7;
					data[datasum].fac=1;
				}
				else
				{
					i+=6;
					data[datasum].fac=0;
				}
			}
			else if(input[i]=='s'&&input[i+1]=='o'&&input[i+2]=='l')
			{
				data[datasum].ope=117;
				i+=5;
			}
			else if(input[i]=='f'&&input[i+1]=='l')
			{
				data[datasum].ope=118;
				i+=5;
			}
			else if(input[i]=='c'&&input[i+1]=='e')
			{
				data[datasum].ope=119;
				i+=4;
			}
			else if(input[i]=='c'&&(input[i+1]=='('||input[i+1]=='['||input[i+1]=='j'))
			{
				data[datasum].ope=120;
				if(input[i+1]=='j')
				{
					data[datasum].fac=1;
					i+=2;
				}
				else
				{
					data[datasum].fac=0;
				    i++;
				}
			}
			else if(input[i]=='r'&&(input[i+1]=='('||input[i+1]=='['))
			{
				data[datasum].ope=122;
				i++;
			}
			else if(input[i]=='u'&&input[i+1]=='p')
			{
				data[datasum].ope=123;
				i+=2;
			}
			else if(input[i]=='d'&&input[i+1]=='o')
			{
				data[datasum].ope=124;
				i+=4;
			}
			else if(solvemode&&input[i]=='q')data[datasum].ope=125;
			else if(input[i]=='c'&&input[i+1]=='o'&&input[i+2]=='d'&&(input[i+5]=='('||input[i+5]=='['))
			{
				data[datasum].ope=126;
				if(input[i+3]=='r'&&input[i+4]=='p')data[datasum].fac=1;
				else if(input[i+3]=='p'&&input[i+4]=='r')data[datasum].fac=2;
				else if(input[i+3]=='r'&&input[i+4]=='s')data[datasum].fac=3;
				else if(input[i+3]=='s'&&input[i+4]=='r')data[datasum].fac=4;
				else if(input[i+3]=='r'&&input[i+4]=='c')data[datasum].fac=5;
				else if(input[i+3]=='c'&&input[i+4]=='r')data[datasum].fac=6;
				else if(input[i+3]=='s'&&input[i+4]=='c')data[datasum].fac=7;
				else if(input[i+3]=='c'&&input[i+4]=='s')data[datasum].fac=8;
				else syntaxerror=1;
				i+=5;
			}
			else if(input[i]=='s'&&input[i+1]=='u'&&input[i+2]=='m'&&(input[i+3]=='X'||input[i+3]=='Y'))
			{
				data[datasum].ope=127;
				if(input[i+3]=='X')
				{
					if(input[i+4]=='Y')
					{
						if(input[i+5]=='['||input[i+5]=='(')
						{
							data[datasum].fac=11;
							i+=5;
						}
						else if(input[i+5]>='1'&&input[i+5]<='9')
						{
							data[datasum].fac=(input[i+5]-'0')*10+1;
							i+=6;
						}
						else syntaxerror=1;
					}
					else if(input[i+4]>='1'&&input[i+4]<='9')
					{
						data[datasum].fac=input[i+4]-'0';
						if(input[i+5]=='['||input[i+5]=='(')i+=5;
						else if(input[i+5]=='Y')
						{
							if(input[i+6]=='['||input[i+6]=='(')
							{
								data[datasum].fac+=10;
								i+=6;
							}
							else if(input[i+6]>='1'&&input[i+6]<='9')
							{
								data[datasum].fac+=10*(input[i+6]-'0');
								i+=7;
							}
							else syntaxerror=1;
						}
						else syntaxerror=1;
					}
					else if(input[i+4]=='['||input[i+4]=='(')
					{
						data[datasum].fac=1;
						i+=4;
					}
					else syntaxerror=1;
				}
				else if(input[i+3]=='Y')
				{
					if(input[i+4]=='['||input[i+4]=='(')
					{
						data[datasum].fac=10;
						i+=4;
					}
					else if(input[i+4]>='1'&&input[i+4]<='9')
					{
						data[datasum].fac=10*(input[i+4]-'0');
						i+=5;
					}
					else syntaxerror=1;
				}
				if(syntaxerror)return;
			}
			else if(input[i]=='f'&&input[i+1]=='c')
			{
				data[datasum].ope=128;
				if(input[i+2]=='('||input[i+2]=='[')
				{
					data[datasum].fac=1;
					i+=2;
				}
				else if(input[i+2]>='1'&&input[i+2]<='9')
				{
					if(input[i+3]=='('||input[i+3]=='[')
					{
						data[datasum].fac=input[i+2]-'0';
						i+=3;
					}
					else if(input[i+3]>='0'&&input[i+3]<='9')
					{
						data[datasum].fac=(input[i+2]-'0')*10+input[i+3]-'0';
						if(data[datasum].fac>43)syntaxerror=1;
						i+=4;
					}
					else syntaxerror=1;
				}
				else if(input[i+2]=='e')
				{
					data[datasum].fac=-1;
					i+=3;
				}
				else if(input[i+2]=='l')
				{
					data[datasum].fac=-2;
					i+=3;
				}
				else if(input[i+2]=='p')
				{
					data[datasum].fac=-3;
					i+=3;
				}
				else if(input[i+2]=='i')
				{
					data[datasum].fac=-4;
					i+=3;
				}
				else syntaxerror=1;
				if(syntaxerror)return;
			}
			else if(input[i]=='A'&&input[i+1]=='[')
			{
				data[datasum].ope=129;
				i++;
			}
			else if(input[i]=='(')data[datasum].ope=11;
			else 
			{
				syntaxerror=1;
				return;
			}
			t=t2=0;
		}
		else if(input[i]==')'||input[i]==']'||input[i]=='}'||(input[i]=='d'&&input[i+1]>='W'&&input[i+1]<='Z'))
		{
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=6;
			if(input[i]==')'||input[i]==']'||input[i]=='}')data[datasum].ope=14;
			else if(input[i]=='d')
			{
				data[datasum].ope=input[i+1]-27>60?input[i+1]-27:64;
				i+=2;
			}
			t=t2=1;
		}
		else if(input[i]==','||input[i]==';'||input[i]==':')
		{
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=7;
			if(input[i]==',')data[datasum].ope=33;
			else if(input[i]==';')data[datasum].ope=99;
			else if(input[i]==':')data[datasum].ope=100;
			t=t2=0;
		}
		else if(input[i]=='e'||(input[i]=='p'&&input[i+1]=='i')||input[i]=='R'||input[i]=='M'||input[i]=='X'||input[i]=='Y'||input[i]=='Z'||(input[i]=='W'&&(input[i+1]<'0'||input[i+1]>'9'))||(input[i]=='A'&&(input[i+1]=='N'||input[i+1]=='T'))||input[i]=='k'||input[i]=='K'||input[i]=='C'||input[i]=='G'||input[i]=='E'||input[i]=='N'||input[i]=='H'||input[i]=='m'||input[i]=='g')
		{
			if(t2==1)
			{
				datasum++;
				data[datasum].id=datasum;
			    data[datasum].cla=2;
				data[datasum].ope=3;
			}
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=1;
			if(input[i]=='e'&&input[i+1]=='0')
			{
				data[datasum].num=E0;
				i++;
			}
			else if(input[i]=='e')data[datasum].num=Ee;
			else if(input[i]=='p')
			{
				data[datasum].num=PI;
				i++;
			}
			else if(input[i]=='R')
			{
				if(input[i+1]=='Y')
				{
					data[datasum].num=RD;
					i++;
				}
				else if(input[i+1]=='O')
				{
					data[datasum].num=Ro;
					i++;
				}
				else if(input[i+1]>='0'&&input[i+1]<='9')
				{
					if(input[i+2]>='0'&&input[i+2]<='9')
					{
						data[datasum].num=randomnumlist[10*(input[i+1]-'0')+input[i+2]-'0'];
						i+=2;
					}
					else
					{
						data[datasum].num=randomnumlist[input[i+1]-'0'];
						i++;
					}
				}
				else data[datasum].num=ER;
			}
			else if(input[i]=='A')
			{
				if(input[i+1]=='N')
				{
					if(memans.fac==0)data[datasum].num=memans.mat[1];
					else
					{
						if(memans.fac==1)data[datasum].fac=1;
						else data[datasum].fac=0;
						data[datasum].cla=16;
						data[datasum].r=memans.r;
						data[datasum].c=memans.c;
						data[datasum].mat=new double[memans.r*memans.c+2];
						for(j=1;j<=memans.r*memans.c;j++)
							data[datasum].mat[j]=memans.mat[j];
					}
				}
				else if(input[i+1]=='T')data[datasum].num=ATM;
				i+=2;
			}
			else if(input[i]=='K')data[datasum].num=Kk;
			else if(input[i]=='C')
			{
				if(input[i+1]=='W'&&input[i+2]=='Y')
				{
					data[datasum].num=CWY;
					i+=2;
				}
				else if(input[i+1]=='O')
				{
					data[datasum].num=Cc;
					i++;
				}
			}
			else if(input[i]=='G')data[datasum].num=Gg;
			else if(input[i]=='E')data[datasum].num=ee;
			else if(input[i]=='N'&&input[i+1]=='A')
			{
				data[datasum].num=NA;
				i++;
			}
			else if(input[i]=='H')data[datasum].num=Hh;
			else if(input[i]=='m')
			{
				if(input[i+1]=='e')data[datasum].num=me;
				else if(input[i+1]=='p')data[datasum].num=mp;
				else if(input[i+1]=='E')data[datasum].num=mE;
				i++;
			}
			else if(input[i]=='g')data[datasum].num=gg;
			else if(input[i]=='M')
			{
				if(input[i+1]>='a'&&input[i+1]<='z')input[i+1]-=32;
				if(input[i+1]>='A'&&input[i+1]<='Z')
				{
					if(mem[input[i+1]-64].fac==0)data[datasum].num=mem[input[i+1]-64].mat[1];
					else
					{
						if(mem[input[i+1]-64].fac==1)data[datasum].fac=1;
						else data[datasum].fac=0;
						data[datasum].cla=16;
						data[datasum].r=mem[input[i+1]-64].r;
						data[datasum].c=mem[input[i+1]-64].c;
						data[datasum].mat=new double[mem[input[i+1]-64].r*mem[input[i+1]-64].c+2];
						for(j=1;j<=mem[input[i+1]-64].r*mem[input[i+1]-64].c;j++)
							data[datasum].mat[j]=mem[input[i+1]-64].mat[j];
					}
				}
				else
				{
					syntaxerror=1;
					return;
				}
				i++;
			}
			else if(input[i]=='X')
			{
				if(input[i+1]<48||input[i+1]>57)data[datasum].cla=9;
				else
				{
					if(input[i+2]>=48&&input[i+2]<=57&&(10*input[i+1]-480+input[i+2]-48)<=solvesum&&(10*input[i+1]-480+input[i+2]-48)>0)
					{
						data[datasum].num=solvelist[10*input[i+1]-480+input[i+2]-48];
						i+=2;
					}
					else if(input[i+1]-48<=solvesum&&input[i+1]>48)
					{
						data[datasum].num=solvelist[input[i+1]-48];
						i++;
					}
					else
					{
						syntaxerror=1;
						return;
					}
				}
			}
			else if(input[i]=='Y')data[datasum].cla=10;
			else if(input[i]=='Z')data[datasum].cla=11;
			else if(input[i]=='W')data[datasum].cla=12;
			else if(input[i]=='k')data[datasum].cla=8;
			t=t2=1;
		}
		else if(input[i]=='-'&&input[i+1]=='>'&&input[i+2]=='M')
		{
			mempos=input[i+3]-64;
			memope=0;
			i+=3;
			if(input[i+1]=='+'||input[i+1]=='-'||input[i+1]=='*'||input[i+1]=='$')
			{
				if(input[i+1]=='+')memope=1;
				else if(input[i+1]=='-')memope=2;
				else if(input[i+1]=='*')memope=3;
				else memope=4;
				i++;
			}
		}
		else if(input[i]=='i'||(input[i]=='W'&&input[i+1]>='0'&&input[i+1]<='9'))
		{
			if(t2==1)
			{
				datasum++;
				data[datasum].id=datasum;
			    data[datasum].cla=2;
				data[datasum].ope=3;
			}
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=16;
			if(input[i]=='i')createcomplex(datasum,0,1);
			else if(input[i]=='W')
			{
				templength=strlen(temp);
				for(j=0;j<=templength;j++)
					temp[j]='\0';
				i++;
				j=i;
				while(input[j]>='0'&&input[j]<='9')
				{
					temp[j-i]=input[j];
					j++;
				}
				i=j-1;
				j=atoi(temp);
				if(j<=0||j>10000)
				{
					syntaxerror=1;
					data[datasum].cla=1;
					return;
				}
				createcomplex(datasum,cos(2*PI/j),sin(2*PI/j));
			}
			t=t2=1;
		}
		else if(input[i]=='T'&&input[i+1]=='O'&&input[i+2]=='D'&&input[i+3]=='A'&&input[i+4]=='Y')
		{
			datasum++;
			data[datasum].id=datasum;
			data[datasum].cla=17;
			SYSTEMTIME sys;
			GetLocalTime(&sys);
			data[datasum].fac=sys.wYear;
			data[datasum].r=sys.wMonth;
			data[datasum].c=sys.wDay;
			i+=4;
			t=t2=1;
		}
		else if((input[i]>=48&&input[i]<=57)||input[i]=='.')
		{
			templength=strlen(temp);
			pointsum=0;
			for(j=0;j<=templength;j++)
				temp[j]='\0';
			j=i;
			while((input[j]>=48&&input[j]<=57)||input[j]=='.')
			{
				temp[j-i]=input[j];
				if(input[j]=='.')pointsum++;
				if(pointsum>=3)
				{
					syntaxerror=1;
					return;
				}
				j++;
			}
			i=j-1;
			datasum++;
			data[datasum].id=datasum;
			if(pointsum<=1)
			{
			    data[datasum].cla=1;
			    data[datasum].num=atof(temp);
			}
			else
			{
				data[datasum].cla=17;
				templength=strlen(temp);
				if(templength>11||templength<5)
				{
					syntaxerror=1;
					return;
				}
				data[datasum].fac=data[datasum].r=data[datasum].c=0;
				ofstream outfile;
				outfile.open("temp.txt",ios::out);
				for(j=0;j<templength;j++)
					outfile<<temp[j];
				outfile.close();
				ifstream infile;
				infile.open("temp.txt",ios::in);
				infile>>data[datasum].fac;
				infile.get();
				infile>>data[datasum].r;
				infile.get();
				infile>>data[datasum].c;
				infile.close();
				remove("temp.txt");
				if(data[datasum].fac>=100000)
				{
					syntaxerror=1;
					return;
				}
			}
			t=t2=1;
		}
		else
		{
			syntaxerror=1;
			return;
		}
		i++;
	}
	addrightbrac();
}
/*----------------------------------------------------------------------------------------------------------*/
void addrightbrac()
{
	int i,tempdatasum,lbs=0,rbs=0;
	for(i=1;i<=datasum;i++)
	{
		if(data[i].cla==5)lbs++;
		if(data[i].cla==6)rbs++;
	}
	tempdatasum=datasum;
	if(lbs>rbs)
	{
	    for(i=tempdatasum+1;i<=tempdatasum+lbs-rbs;i++)
		{
		    datasum++;
		    data[datasum].id=datasum;
		    data[datasum].cla=6;
		    data[datasum].ope=14;
		}
	}
	else if(lbs<rbs)
	{
		for(i=tempdatasum;i>tempdatasum-rbs+lbs;i--)
		{
			data[datasum].id=0;
			datasum--;
		}
	}
}
/*----------------------------------------------------------------------------------------------------------*/
int searchbrac()
{
	int i,r=11;
	ismatrix=0;
	for(i=1;i<=datasum;i++)
	{
		if(data[i].id>0&&data[i].cla==5)
		{
			r=data[i].ope;
			if(r==37||r==38||r==39||r==96)
			{
				calcexpression(i,r);
				r=11;
			}
			leftbrac=i;
		}
		if(data[i].id>0&&data[i].cla==6)
		{
			rightbrac=i;
			break;
		}
	}
	data[leftbrac].id=data[rightbrac].id=0;
	if(leftbrac>0&&rightbrac>0)dataleft-=2;
	if(r==11||r==13||(r>=16&&r<=18)||r==21||r==34||r==36||r==55||r==56||r==113||r==20||r==28||(r>=25&&r<=27)||r==118||r==119||r==125||r==72)
	{
		for(i=leftbrac+1;i<rightbrac;i++)
			if(data[i].id>0&&(data[i].cla==7||data[i].cla==16))
			{
				ismatrix=1;
				break;
			}
	}
	return r;
}
/*----------------------------------------------------------------------------------------------------------*/
void calcexpression(const int x,const int z)
{
	int rightbracpos,comma1pos,comma2pos,i,j,k,dts=datasum,dtl=dataleft,leftbracsum=0,rightbracsum=0,wtime,tempme,tempme2=0,dvar=1,idivideparttemp,tempi;
	double result,ia,ib,h,tempt,resultim=0;
	bool find,mat=0,cmplx=1,matme=0,calcim=0;
	__int64 startk,endk,tempk;
	countfeverexpressiondata *tempcfr=new countfeverexpressiondata[dts+2];
	if(z==39)idivideorder++;
	if(idivideorder==1)
	{
		if(data[x].fac!=-1)idivideparttemp=idividepart[data[x].fac];
		else idivideparttemp=40000;
	}
	else
	{
		if(data[x].fac>=idivideorder)idivideparttemp=idividepart[data[x].fac];
		else idivideparttemp=idividepart[idivideorder];
	}
	if(funmode&&idivideparttemp>2000)idivideparttemp=2000;
	if((solvemode||plotting)&&idivideparttemp>200)idivideparttemp=200;
	if(z==38)result=1;
	else result=0;
	for(i=1;i<=dts;i++)
	{
		if(data[i].cla!=16)tempcfr[i]=data[i];
		else
		{
			tempcfr[i].r=data[i].r;
			tempcfr[i].c=data[i].c;
			tempcfr[i].fac=data[i].fac;
			tempcfr[i].id=data[i].id;
			tempcfr[i].num=data[i].num;
			tempcfr[i].cla=16;
			tempcfr[i].mat=new double[tempcfr[i].r*tempcfr[i].c+2];
			for(j=1;j<=tempcfr[i].r*tempcfr[i].c;j++)
				tempcfr[i].mat[j]=data[i].mat[j];
		}
	}
	for(i=x+1; ;i++)
	{
		if(tempcfr[i].cla==5)leftbracsum++;
		if(tempcfr[i].cla==6)rightbracsum++;
		if(rightbracsum==leftbracsum&&tempcfr[i].cla==7)break;
		if(i>dts)
		{
			syntaxerror=1;
			break;
		}
	}
	comma1pos=i;
	datasum=dataleft=i-x-1;
	for(i=x+1;i<comma1pos;i++)
	{
		if(tempcfr[i].cla!=16)data[i-x]=tempcfr[i];
		else
		{
			data[i-x].cla=16;
			data[i-x].r=tempcfr[i].r;
			data[i-x].c=tempcfr[i].c;
			data[i-x].fac=tempcfr[i].fac;
			data[i-x].id=tempcfr[i].id;
			data[i-x].mat=new double[tempcfr[i].r*tempcfr[i].c+2];
			for(j=1;j<=tempcfr[i].r*tempcfr[i].c;j++)
			    data[i-x].mat[j]=tempcfr[i].mat[j];
		}
	}
	dtl-=(datasum-1);
	wtime=maincalc();
	find=0;
	if(syntaxerror==0&&wtime<400)
		for(i=1;i<=datasum;i++)
		{
			if(data[i].id>0&&data[i].cla==16)
			{
				if(data[i].r==1&&data[i].c==1)oneomatrix(i);
				else if(data[i].r==1&&data[i].c==2&&fabs(data[i].mat[2])<1e-15)oneomatrix(i);
				else
				{
					delete []data[i].mat;
					matherr(i);
				}
			}
			if(data[i].id>0&&_isnan(data[i].num)==0)
			{
				if(z==39)
				{
					ia=data[i].num;
					startk=0;
				}
				else if(z==96)
				{
					startk=0;
					endk=4;
					if(fabs(data[i].num)<1e-4)
					{
						ia=data[i].num-(2e-11);
						ib=data[i].num+(2e-11);
						h=1e-11;
					}
					else
					{
						ia=data[i].num-fabs(data[i].num)/(1e+7);
						ib=data[i].num+fabs(data[i].num)/(1e+7);
						h=fabs(data[i].num)/(2e+7);
					}
				}
				else startk=__int64(data[i].num);
				find=1;
				break;
			}
		}
	if(find==0)syntaxerror=1;
    if(z>=37&&z<=39)
	{
	    for(i=comma1pos+1; ;i++)
		{
		    if(tempcfr[i].cla==5)leftbracsum++;
		    if(tempcfr[i].cla==6)rightbracsum++;
		    if(rightbracsum==leftbracsum&&tempcfr[i].cla==7)break;
			if(i>dts)
			{
			    syntaxerror=1;
			    break;
			}
		}
	    comma2pos=i;
	    datasum=dataleft=i-comma1pos-1;
	    for(i=comma1pos+1;i<comma2pos;i++)
		{
		    if(tempcfr[i].cla!=16)data[i-comma1pos]=tempcfr[i];
	    	else
			{
			    data[i-comma1pos].cla=16;
			    data[i-comma1pos].r=tempcfr[i].r;
			    data[i-comma1pos].c=tempcfr[i].c;
			    data[i-comma1pos].fac=tempcfr[i].fac;
			    data[i-comma1pos].id=tempcfr[i].id;
			    data[i-comma1pos].mat=new double[tempcfr[i].r*tempcfr[i].c+2];
			    for(j=1;j<=tempcfr[i].r*tempcfr[i].c;j++)
			        data[i-comma1pos].mat[j]=tempcfr[i].mat[j];
			}
		}
	    dtl-=(datasum-1);
	    wtime=maincalc();
	    find=0;
	    if(syntaxerror==0&&wtime<400)
		    for(i=1;i<=datasum;i++)
			{
				if(data[i].id>0&&data[i].cla==16)
				{
				    if(data[i].r==1&&data[i].c==1)oneomatrix(i);
				    else if(data[i].r==1&&data[i].c==2&&fabs(data[i].mat[2])<1e-15)oneomatrix(i);
				    else
					{
					    delete []data[i].mat;
					    matherr(i);
					}
				}
			    if(data[i].id>0&&_isnan(data[i].num)==0)
				{
				    if(z==39)
					{
					    ib=data[i].num;
					    endk=idivideparttemp;
					}
				    else endk=__int64(data[i].num);
				    find=1;
				    break;
				}
			}
	    if(find==0)syntaxerror=1;
	}
	else if(z==96)comma2pos=comma1pos;
	for(i=comma2pos+1; ;i++)
	{
		if(tempcfr[i].cla==5)leftbracsum++;
		if(tempcfr[i].cla==6)rightbracsum++;
		if(rightbracsum-leftbracsum==1)break;
		if(i>dts)
		{
		    syntaxerror=1;
		    break;
		}
	}
	rightbracpos=i;
	if(tempcfr[i].ope>=62)dvar=tempcfr[i].ope-60;
	datasum=dataleft=i-comma2pos-1;
	dtl-=(datasum-1);
	tempme=matherror;
	for(tempk=startk;tempk<=endk;tempk++)
	{
		for(i=comma2pos+1;i<rightbracpos;i++)
		{
		    if(tempcfr[i].cla!=16)data[i-comma2pos]=tempcfr[i];
		    else
			{
			    data[i-comma2pos].cla=16;
			    data[i-comma2pos].r=tempcfr[i].r;
			    data[i-comma2pos].c=tempcfr[i].c;
			    data[i-comma2pos].fac=tempcfr[i].fac;
			    data[i-comma2pos].id=tempcfr[i].id;
			    data[i-comma2pos].mat=new double[tempcfr[i].r*tempcfr[i].c+2];
			    for(j=1;j<=tempcfr[i].r*tempcfr[i].c;j++)
				    data[i-comma2pos].mat[j]=tempcfr[i].mat[j];
			}
			if(z<=38&&data[i-comma2pos].cla==8)
			{
				data[i-comma2pos].cla=1;
				data[i-comma2pos].num=tempk;
			}
			else if((z==39||z==96)&&data[i-comma2pos].cla==dvar+8)
			{
				data[i-comma2pos].cla=1;
				if(z==39)data[i-comma2pos].num=ia+(ib-ia)*tempk/idivideparttemp;
				else data[i-comma2pos].num=ia+(ib-ia)*tempk/endk;
			}
		}
		dataleft=rightbracpos-comma2pos-1;
		wtime=maincalc();
		if(syntaxerror==0&&wtime<400)
			for(i=1;i<=datasum;i++)
				if(data[i].id>0&&((data[i].cla!=16&&_isnan(data[i].num)==0)||(data[i].cla==16)))
				{
					if(z==37)
					{
						if(data[i].cla!=16||(data[i].cla==16&&data[i].r==1&&data[i].c==2))
						{
							if(!mat)
							{
								if(data[i].cla!=16)result+=data[i].num;
								else
								{
									if(!calcim)calcim=1;
									result+=data[i].mat[1];
									resultim+=data[i].mat[2];
									delete []data[i].mat;
								}
							}
							else matme=1;
						}
						else
						{
							if(!mat)
							{
								mat=1;
								if(result!=0)matme=1;
								if(data[i].fac==0)cmplx=0;
								data[0].r=data[i].r;
								data[0].c=data[i].c;
								data[0].fac=data[i].fac;
								data[0].mat=new double[data[0].r*data[0].c+2];
								for(j=0;j<=data[0].r*data[0].c;j++)
									data[0].mat[j]=0;
							}
							if(data[i].fac==0)cmplx=0;
							countmatrix(0,i,1);
							if(matherror)matme=1;
						}
					}
					else if(z==38)
					{
						if(data[i].cla!=16||(data[i].cla==16&&data[i].r==1&&data[i].c==2))
						{
							if(!mat)
							{
								if(data[i].cla!=16)
								{
									result*=data[i].num;
									resultim*=data[i].num;
								}
								else
								{
									if(!calcim)calcim=1;
									double tempre=result*data[i].mat[1]-resultim*data[i].mat[2],tempim=result*data[i].mat[2]+resultim*data[i].mat[1];
									result=tempre;
									resultim=tempim;
									delete []data[i].mat;
								}
							}
							else matme=1;
						}
						else
						{
							if(!mat)
							{
								mat=1;
								if(result!=1)matme=1;
								if(data[i].fac==0)cmplx=0;
								data[0].r=data[0].c=data[i].r;
								data[0].fac=data[i].fac;
								data[0].mat=new double[data[0].r*data[0].c+2];
								for(j=0;j<=data[0].r*data[0].c;j++)
									data[0].mat[j]=0;
								for(j=1;j<=data[0].r*data[0].c;j+=data[0].r+1)
									data[0].mat[j]=1;
							}
							if(data[i].fac==0)cmplx=0;
							countmatrix(0,i,3);
							if(matherror)matme=1;
						}
					}
					else if(z==39)
					{
						if(tempk%4==0&&(tempk==startk||tempk==endk))tempi=7;
						else if(tempk%4==0)tempi=14;
					    else if(tempk%4==1||tempk%4==3)tempi=32;
					    else tempi=12;
						if(data[i].cla!=16||(data[i].cla==16&&data[i].r==1&&data[i].c==2))
						{
							if(!mat)
							{
								if(data[i].cla!=16)result+=tempi*data[i].num;
								else
								{
									if(!calcim)calcim=1;
									result+=tempi*data[i].mat[1];
									resultim+=tempi*data[i].mat[2];
									delete []data[i].mat;
								}
							}
							else matme=1;
						}
						else
						{
							if(!mat)
							{
								mat=1;
								if(result!=0)matme=1;
								if(data[i].fac==0)cmplx=0;
								data[0].r=data[i].r;
								data[0].c=data[i].c;
								data[0].fac=data[i].fac;
								data[0].mat=new double[data[0].r*data[0].c+2];
								for(j=0;j<=data[0].r*data[0].c;j++)
									data[0].mat[j]=0;
							}
							if(data[i].fac==0)cmplx=0;
							for(k=1;k<=data[i].r*data[i].c;k++)
								data[i].mat[k]*=tempi;
							countmatrix(0,i,1);
							if(matherror)matme=1;
						}
					}
					else if(z==96)
					{
						if(tempk==0)tempi=1;
						else if(tempk==1)tempi=-8;
						else if(tempk==2)tempi=0;
						else if(tempk==3)tempi=8;
						else if(tempk==4)tempi=-1;
						if(data[i].cla!=16||(data[i].cla==16&&data[i].r==1&&data[i].c==2))
						{
							if(!mat)
							{
								if(data[i].cla!=16)result+=tempi*data[i].num;
								else
								{
									if(!calcim)calcim=1;
									result+=tempi*data[i].mat[1];
									resultim+=tempi*data[i].mat[2];
									delete []data[i].mat;
								}
							}
							else matme=1;
						}
						else
						{
							if(!mat)
							{
								mat=1;
								if(result!=0)matme=1;
								if(data[i].fac==0)cmplx=0;
								data[0].r=data[i].r;
								data[0].c=data[i].c;
								data[0].fac=data[i].fac;
								data[0].mat=new double[data[0].r*data[0].c+2];
								for(j=0;j<=data[0].r*data[0].c;j++)
									data[0].mat[j]=0;
							}
							if(data[i].fac==0)cmplx=0;
							for(k=1;k<=data[i].r*data[i].c;k++)
								data[i].mat[k]*=tempi;
							countmatrix(0,i,1);
							if(matherror)matme=1;
						}
					}
					tempme2=1;
					break;
				}
	    if(z==39&&tempme==0)matherror=0;
	}
	if((z==39&&tempme2==0)||(mat&&matme))matherror=1;
	if(z>=37&&z<=39)dataleft=dtl-4;
	else dataleft=dtl-2;
	datasum=dts;
	for(i=1;i<=datasum;i++)
	{
		if(tempcfr[i].cla!=16)data[i]=tempcfr[i];
		else
		{
			data[i].cla=16;
			data[i].id=tempcfr[i].id;
			data[i].r=tempcfr[i].r;
			data[i].c=tempcfr[i].c;
			data[i].fac=tempcfr[i].fac;
			data[i].mat=new double[data[i].r*data[i].c+2];
			for(j=1;j<=data[i].r*data[i].c;j++)
				data[i].mat[j]=tempcfr[i].mat[j];
		}
		if(i>x&&i<rightbracpos)data[i].id=0;
	}
	if(z==39&&!mat)
	{
		result*=(2*(ib-ia)/idivideparttemp/45);
		if(calcim)resultim*=(2*(ib-ia)/idivideparttemp/45);
		if(fabs(result)<1e-15)result=0;
	}
	else if(z==39)
	{
		tempt=2*(ib-ia)/idivideparttemp/45;
		for(k=1;k<=data[0].r*data[0].c;k++)
		{
			data[0].mat[k]*=tempt;
			if(fabs(data[0].mat[k])<1e-15)result=0;
		}
	}
	else if(z==96&&!mat)
	{
		result/=(12*h);
		if(calcim)resultim/=(12*h);
	}
	else if(z==96)
	{
		tempt=12*h;
		for(k=1;k<=data[0].r*data[0].c;k++)
			data[0].mat[k]/=tempt;
	}
	data[x+1].id=x+1;
	if(!mat)
	{
		if(!calcim||(calcim&&fabs(resultim)<1e-12))
		{
	        data[x+1].cla=1;
	        data[x+1].num=result;
		}
		else createcomplex(x+1,result,resultim);
	}
	else
	{
		data[x+1].cla=16;
		data[x+1].r=data[0].r;
		data[x+1].c=data[0].c;
		data[x+1].fac=cmplx;
		data[x+1].mat=new double[data[x+1].r*data[x+1].c+2];
		for(i=1;i<=data[x+1].r*data[x+1].c;i++)
			data[x+1].mat[i]=data[0].mat[i];
		delete []data[0].mat;
	}
	if(z==39)idivideorder--;
	for(i=1;i<=dts;i++)
		if(tempcfr[i].cla==16)delete []tempcfr[i].mat;
	delete []tempcfr;
}
/*----------------------------------------------------------------------------------------------------------*/
void calcvector(const int x,const int y)
{
	int z=data[x].ope,firstnum=0,i,j,k,commasum=0,lastpos=x,numsum=0,sp=data[x].fac;
	double optinum,tempnum,tempnum2,number[2001]={0},sumofset=0,sum2ofset=0,sumofreset=0,proofset=1,average,numsumf,a,b,c,d,r,arg;
	bool find;
	data[y].id=y;
	for(i=x+1;i<=y;i++)
		if(data[i].id>0&&(data[i].cla==7||data[i].cla==6))
		{
			data[i].id=0;
			calculate(lastpos,i,20);
			find=0;
			for(j=lastpos+1;j<i;j++)
			{
				if(data[j].id>0&&data[j].cla==1)
				{
					data[j].id=0;
					tempnum=data[j].num;
					numsum++;
					if(numsum>=2000)
					{
						matherror=1;
						return;
					}
					number[numsum]=tempnum;
					sumofset+=tempnum;
					sum2ofset+=tempnum*tempnum;
					proofset*=tempnum;
					sumofreset+=1/tempnum;
					find=1;
					break;
				}
				else if(data[j].id>0&&data[j].cla==16)
				{
					data[j].id=0;
					for(k=1;k<=data[j].c*data[j].r;k++)
					{
						numsum++;
						tempnum=data[j].mat[k];
						if(numsum>=2000)
						{
							matherror=1;
							delete []data[j].mat;
							return;
						}
						number[numsum]=tempnum;
					    sumofset+=tempnum;
					    sum2ofset+=tempnum*tempnum;
				    	proofset*=tempnum;
					    sumofreset+=1/tempnum;
						if(firstnum==0)
						{
				            firstnum=1;
				            optinum=tempnum;
						}
						else
						{
				            if(z==43&&tempnum>optinum)optinum=tempnum;
				            else if(z==44&&tempnum<optinum)optinum=tempnum;
				            else if(z==31)optinum=gcd(optinum,tempnum,&matherror);
				            else if(z==32)optinum=lcm(optinum,tempnum,&matherror);
						}
					}
					find=1;
					delete []data[j].mat;
					break;
				}
			}
			if(find==0)syntaxerror=1;
			if(data[i].cla==7)commasum++;
			lastpos=i;
			if(firstnum==0)
			{
				firstnum=1;
				optinum=tempnum;
			}
			else
			{
				if(z==43&&tempnum>optinum)optinum=tempnum;
				else if(z==44&&tempnum<optinum)optinum=tempnum;
				else if(z==31)optinum=gcd(optinum,tempnum,&matherror);
				else if(z==32)optinum=lcm(optinum,tempnum,&matherror);
			}
		}
	dataleft-=(2*commasum);
	if(z==73)
	{
		if(numsum!=4&&numsum!=5&&numsum!=7)
		{
			matherror=1;
			optinum=0;
		}
		else if(numsum==4)
		{
			if(number[1]==0)matherror=1;
			number[4]=number[4]>=0?1:-1;
			if(number[2]*number[2]-4*number[1]*number[3]>=0)optinum=(-number[2]+number[4]*sqrt(number[2]*number[2]-4*number[1]*number[3]))/(2*number[1]);
			else
			{
				createcomplex(x,-number[2]/(2*number[1]),(number[4]*sqrt(4*number[1]*number[3]-number[2]*number[2]))/(2*number[1]));
				return;
			}
		}
		else if(numsum==5)
		{
			if(number[1]==0)matherror=1;
			else if(number[1]<0)
			{
				for(i=1;i<=5;i++)
					number[i]*=-1;
			}
			double delta,A=number[2]*number[2]-3*number[1]*number[3],B=number[2]*number[3]-9*number[1]*number[4],C=number[3]*number[3]-3*number[2]*number[4],K,T;
			if(number[5]!=0)number[5]=number[5]>0?1:-1;
			K=B/A;
			a=number[1];
			b=number[2];
			c=number[3];
			d=number[4];
			T=(2*A*b-3*a*B)/(2*A*sqrt(A));
			delta=B*B-4*A*C;
			if(b==0&&c==0&&d==0)optinum=0;
			else if(A==0&&B==0)optinum=-c/b;
			else if(delta==0)
			{
				if(number[5]==0)optinum=-b/a+K;
				else optinum=-K/2;
			}
			else if(delta<0)
			{
				double tempdb=acos(T)/3;
				if(number[5]==0)optinum=(-b-2*sqrt(A)*cos(tempdb))/(3*a);
				else optinum=(-b+sqrt(A)*(cos(tempdb)+number[5]*sqrt(3)*sin(tempdb)))/(3*a);
			}
			else
			{
				double Y1=A*b+3*a*(-B+sqrt(delta))/2,Y2=A*b+3*a*(-B-sqrt(delta))/2;
				if(number[5]==0)optinum=(-b-__pow(Y1,1.0/3,&matherror)-__pow(Y2,1.0/3,&matherror))/(3*a);
				else
				{
					createcomplex(x,(-2*b+__pow(Y1,1.0/3,&matherror)+__pow(Y2,1.0/3,&matherror))/(6*a),number[5]*sqrt(3)*(__pow(Y1,1.0/3,&matherror)-__pow(Y2,1.0/3,&matherror))/(6*a));
					return;
				}
			}
		}
		else if(numsum==7)
		{
			if(number[1]==0&&number[2]==0)matherror=1;
			number[7]=number[7]>=0?1:-1;
			a=number[3]*number[3]-number[4]*number[4]-4*number[1]*number[5]+4*number[2]*number[6];
			b=2*number[3]*number[4]-4*number[2]*number[5]-4*number[1]*number[6];
			c=sqrt(sqrt(a*a+b*b));
			d=planeangle(a,b)/2;
			if(d>PI/2)d-=PI;
			a=-number[3]+number[7]*c*cos(d);
			b=-number[4]+number[7]*c*sin(d);
			c=2*number[1];
			d=2*number[2];
			createcomplex(x,(a*c+b*d)/(c*c+d*d),(b*c-a*d)/(c*c+d*d));
			return;
		}
	}
	else if(z==74)
	{
		if(numsum==1)
		{
			matherror=1;
			optinum=0;
		}
		else if(numsum%2==1)
		{
			optinum=0;
			for(i=numsum/2+1;i<numsum;i++)
				optinum+=number[i]*number[i];
			if(optinum==0)optinum=matherror=1;
			tempnum=0;
			for(i=1;i<=numsum/2;i++)
				tempnum+=number[i]*number[i+numsum/2];
			optinum=fabs((tempnum+number[numsum])/sqrt(optinum));
		}
		else
		{
			optinum=0;
			for(i=1;i<=numsum/2;i++)
				optinum+=(number[i]-number[i+numsum/2])*(number[i]-number[i+numsum/2]);
			optinum=sqrt(optinum);
		}
	}
	else if(z==84)
	{
		if(sp==0)optinum=sumofset;
		else
		{
			if(numsum%2==1)matherror=1;
			optinum=0;
			for(i=1;i<=numsum;i+=2)
			{
				if(number[i+1]<0)matherror=1;
				else optinum+=number[i]*number[i+1];
			}
		}
	}
	else if(z==85)
	{
		if(sp<10)
		{
			if(sp==1)optinum=sumofset;
			else if(sp==2)optinum=sum2ofset;
			else
			{
				optinum=0;
				for(i=1;i<=numsum;i++)
				{
					tempnum=1;
					for(j=1;j<=sp;j++)
						tempnum*=number[i];
					optinum+=tempnum;
				}
			}
		}
		else
		{
			if(numsum%2==1)matherror=1;
			optinum=0;
			for(i=1;i<=numsum;i+=2)
			{
				if(number[i+1]<=0)matherror=1;
				else
				{
					tempnum=number[i+1];
					for(j=1;j<=sp-10;j++)
						tempnum*=number[i];
					optinum+=tempnum;
				}
			}
		}
	}
	else if(z==86||z==87||z==88||z==90||z==91||z==93)
	{
		numsumf=numsum;
		if(z==87&&numsum==1)matherror=1;
		if(sp==0)average=sumofset/numsum;
		else
		{
			if(numsum%2==1)matherror=1;
			sumofset=sum2ofset=sumofreset=numsumf=0;
			for(i=1;i<=numsum;i+=2)
			{
				if(sp==3&&i==1)i=2;
				if(sp==2||sp==3)number[i+1]=1;
				if(number[i+1]<0)matherror=1;
				else
				{
					sumofset+=number[i]*number[i+1];
					sum2ofset+=number[i]*number[i]*number[i+1];
					numsumf+=number[i+1];
					sumofreset+=number[i+1]/number[i];
				}
			}
			average=sumofset/numsumf;
		}
		if(z==86)optinum=sqrt((sum2ofset-2*average*sumofset+numsumf*average*average)/numsumf);
		else if(z==87)optinum=sqrt((sum2ofset-2*average*sumofset+numsumf*average*average)/(numsumf-1));
		else if(z==88)optinum=average;
		else if(z==90)optinum=numsumf;
		else if(z==91)optinum=sqrt(sum2ofset/numsumf);
		else if(z==93)optinum=numsumf/sumofreset;
	}
	else if(z==89&&sp==0)
	{
		for(i=1;i<=numsum;i++)
			for(j=i+1;j<=numsum;j++)
				if(number[i]-number[j]<0)
				{
					tempnum=number[i];
					number[i]=number[j];
					number[j]=tempnum;
				}
		if(numsum%2==1)optinum=number[(numsum+1)/2];
		else optinum=(number[numsum/2]+number[numsum/2+1])/2;
	}
	else if(z==89&&sp==1)
	{
		if(numsum%2==1)matherror=1;
		for(i=1;i<=numsum;i+=2)
			for(j=i+2;j<=numsum;j+=2)
				if(number[i]-number[j]<0)
				{
					tempnum=number[i];
					number[i]=number[j];
					number[j]=tempnum;
					tempnum=number[i+1];
					number[i+1]=number[j+1];
					number[j+1]=tempnum;
				}
		numsumf=0;
		for(i=2;i<=numsum;i+=2)
		{
			if(number[i]<0)matherror=1;
			numsumf+=number[i];
		}
		tempnum=0;
		for(i=2;i<=numsum;i+=2)
		{
			tempnum+=number[i];
			if(tempnum>numsumf/2)break;
		}
		optinum=number[i-1];
	}
	else if(z==92&&sp==0)
	{
		if(numsum%2==0&&proofset<0)matherror=1;
		else optinum=__pow(proofset,1.0/numsum,&matherror);
	}
	else if((z==92||z==94)&&sp==1)
	{
		if(numsum%2!=0)matherror=1;
		numsumf=0;
		proofset=1;
		for(i=1;i<=numsum;i+=2)
		{
			if(number[i+1]<0)matherror=1;
			else proofset*=__pow(number[i],number[i+1],&matherror);
			numsumf+=number[i+1];
		}
		if(z==92)optinum=__pow(proofset,1.0/numsumf,&matherror);
		else optinum=proofset;
	}
	else if(z==94&&sp==0)optinum=proofset;
	else if(z==97||z==106)
	{
		if((z==106&&!isint(sqrt(numsum)))||(z==97&&!isint(sqrt(numsum))&&(!isint(sqrt(numsum/2))||numsum%2!=0)))
		{
			matherror=1;
			optinum=0;
		}
		else if(z==97&&isint(sqrt(numsum)))optinum=det(int(gauss(sqrt(numsum))),number);
		else if(z==97)
		{
			detcm(x,int(gauss(sqrt(numsum/2))),number);
			return;
		}
		else if(z==106)optinum=tr(int(gauss(sqrt(numsum))),number);
	}
	else if(z==13)optinum=sqrt(sum2ofset);
	else if(z==36)
	{
		if(numsum!=1&&numsum!=3&&numsum!=2)syntaxerror=1;
		else if(numsum==1)optinum=ran(number[1]);
		else if(int(number[1])>0&&int(number[2])>0&&int(number[1])*int(number[2])<=2000&&(numsum==3||numsum==2))
		{
			data[x].id=x;
			data[x].cla=16;
			data[x].fac=0;
			data[x].r=int(number[1]);
			data[x].c=int(number[2]);
			tempnum=numsum==3?number[3]:0;
			data[x].mat=new double[data[x].r*data[x].c+2];
			for(i=1;i<=data[x].r*data[x].c;i++)
				data[x].mat[i]=ran(tempnum);
			return;
		}
		else matherror=1;
	}
	else if(z==113||z==18||z==16||z==17||z==21||z==55||z==56)
	{
		if(numsum%2==1)
		{
			numsum++;
			number[numsum]=0;
		}
		if(numsum==2&&z==113)
		{
			optinum=planeangle(number[1],number[2]);
			if(optinum>PI)optinum-=2*PI;
			if(degreemode)optinum*=(180/PI);
		}
		else if(numsum==2&&number[2]==0)optinum=count(number[1],NULL,z,&matherror,degreemode);
		else if(numsum==2&&number[2]!=0)
		{
			if(degreemode)
			{
				number[1]*=(PI/180);
				number[2]*=(PI/180);
			}
			a=sin(number[1])*cosh(number[2]);
			b=cos(number[1])*sinh(number[2]);
			c=cos(number[1])*cosh(number[2]);
			d=-sin(number[1])*sinh(number[2]);
			createcomplex(x,1,1);
			if(z==16)
			{
				data[x].mat[1]=a;
				data[x].mat[2]=b;
			}
			else if(z==17)
			{
				data[x].mat[1]=c;
				data[x].mat[2]=d;
			}
			else if(z==18)
			{
				data[x].mat[1]=(a*c+b*d)/(c*c+d*d);
				data[x].mat[2]=(b*c-a*d)/(c*c+d*d);
			}
			else if(z==21)
			{
				data[x].mat[1]=(a*c+b*d)/(a*a+b*b);
				data[x].mat[2]=(a*d-b*c)/(a*a+b*b);
			}
			else if(z==55)
			{
				data[x].mat[1]=c/(c*c+d*d);
				data[x].mat[2]=-d/(c*c+d*d);
			}
			else if(z==56)
			{
				data[x].mat[1]=a/(a*a+b*b);
				data[x].mat[2]=-b/(a*a+b*b);
			}
			return;
		}
		else
		{
			optinum=0;
			for(i=1;i<=numsum/2;i++)
				optinum+=number[i]*number[i+numsum/2];
			tempnum=0;
			for(i=1;i<=numsum/2;i++)
				tempnum+=number[i]*number[i];
			tempnum=sqrt(tempnum);
			tempnum2=0;
			for(i=numsum/2+1;i<=numsum;i++)
				tempnum2+=number[i]*number[i];
			tempnum2=sqrt(tempnum2);
			optinum/=(tempnum*tempnum2);
			optinum=count(optinum,NULL,23,&matherror,degreemode);
			if(z!=113)optinum=count(optinum,NULL,z,&matherror,degreemode);
		}
	}
	else if(z==15)
	{
		if(numsum>2)matherror=1;
		if(numsum==1)
		{
			if(number[1]>=0)optinum=sqrt(number[1]);
			else
			{
				createcomplex(x,0,sqrt(-number[1]));
				return;
			}
		}
		else if(numsum==2)
		{
			tempnum=sqrt(sqrt(number[1]*number[1]+number[2]*number[2]));
			tempnum2=planeangle(number[1],number[2])/2;
			if(tempnum2>PI/2)tempnum2-=PI;
			createcomplex(x,tempnum*cos(tempnum2),tempnum*sin(tempnum2));
			return;
		}
	}
	else if(z==12||z==19)
	{
		if(numsum>2)matherror=1;
		else if(numsum==1)
		{
			if(number[1]==0)matherror=1;
			else if(number[1]>0)
			{
				if(z==12)optinum=log(number[1]);
				else optinum=log(number[1])/log(10);
			}
			else
			{
				createcomplex(x,log(-number[1]),PI);
				if(z==19)
				{
					data[x].mat[1]/=log(10);
					data[x].mat[2]/=log(10);
				}
				return;
			}
		}
		else if(numsum==2)
		{
			if(number[1]==0&&number[2]==0)matherror=1;
			else
			{
				createcomplex(x,log(number[1]*number[1]+number[2]*number[2])/2,planeangle(number[1],number[2]));
				if(data[x].mat[2]>PI)data[x].mat[2]-=2*PI;
				if(z==19)
				{
					data[x].mat[1]/=log(10);
					data[x].mat[2]/=log(10);
				}
				return;
			}
		}
	}
	else if(z==20||(z>=25&&z<=27))
	{
		if(numsum==1)optinum=count(number[1],NULL,z,&matherror,degreemode);
		else if(numsum!=2)matherror=1;
		else
		{
			createcomplex(x,1,1);
			if(z==20)
			{
			    data[x].mat[1]=exp(number[1])*cos(number[2]);
			    data[x].mat[2]=exp(number[1])*sin(number[2]);
			}
			else
			{
				a=(exp(number[1])*cos(number[2])-exp(-number[1])*cos(number[2]))/2;
				b=(exp(number[1])*sin(number[2])+exp(-number[1])*sin(number[2]))/2;
				if(z==25)
				{
					data[x].mat[1]=a;
					data[x].mat[2]=b;
				}
				else
				{
					c=(exp(number[1])*cos(number[2])+exp(-number[1])*cos(number[2]))/2;
				    d=(exp(number[1])*sin(number[2])-exp(-number[1])*sin(number[2]))/2;
					if(z==26)
					{
						data[x].mat[1]=c;
						data[x].mat[2]=d;
					}
					else
					{
						data[x].mat[1]=(a*c+b*d)/(c*c+d*d);
						data[x].mat[2]=(b*c-a*d)/(c*c+d*d);
					}
				}
			}
			return;
		}
	}
	else if(z==28)
	{
		if(numsum==1||numsum>4)matherror=1;
		else if(numsum==2&&number[1]>0&&number[2]>0)optinum=count(number[1],number[2],z,&matherror,degreemode);
		else
		{
			if(numsum==2)
			{
				a=number[1];
				c=number[2];
				b=d=0;
			}
			else if(numsum==3)
			{
				a=number[1];
				b=0;
				c=number[2];
				d=number[3];
			}
			else
			{
				a=number[1];
				b=number[2];
				c=number[3];
				d=number[4];
			}
			if((a==0||a==1&&b==0)||(c==0&&d==0))matherror=1;
			else
			{
				createcomplex(x,log(c*c+d*d)/2,planeangle(c,d));
				if(data[x].mat[2]>PI)data[x].mat[2]-=2*PI;
				c=data[x].mat[1];
				d=data[x].mat[2];
				data[x].mat[1]=log(a*a+b*b)/2;
				data[x].mat[2]=planeangle(a,b);
				if(data[x].mat[2]>PI)data[x].mat[2]-=2*PI;
				a=data[x].mat[1];
				b=data[x].mat[2];
				data[x].mat[1]=(a*c+b*d)/(a*a+b*b);
				data[x].mat[2]=(a*d-b*c)/(a*a+b*b);
				return;
			}
		}
	}
	else if(z==114||z==115)
	{
		if(numsum==1)
		{
			if(z==114)optinum=number[1];
			else optinum=0;
		}
		else if(numsum==2)
		{
			if(z==114)optinum=number[1];
			else optinum=number[2];
		}
		else matherror=1;
	}
	else if((z>=22&&z<=24)||z==83)
	{
		if(numsum>2)matherror=1;
		else if(numsum==2)
		{
			c=number[1];
			d=number[2];
		}
		else
		{
			c=number[1];
			d=0;
		}
		if(numsum<=2)
		{
			if(z==22)
			{
				r=sqrt(sqrt((1-c*c+d*d)*(1-c*c+d*d)+4*c*c*d*d));
				arg=planeangle(1-c*c+d*d,-2*c*d)/2;
				if(arg>PI/2)arg-=PI;
				a=r*cos(arg)-d;
				b=r*sin(arg)+c;
			}
			else if(z==23)
			{
				r=sqrt(sqrt((1-c*c+d*d)*(1-c*c+d*d)+4*c*c*d*d));
				arg=planeangle(c*c-d*d-1,2*c*d)/2;
				if(arg>PI/2)arg-=PI;
				a=r*cos(arg)+c;
				b=r*sin(arg)+d;
			}
			else if(z==24||z==83)
			{
				if(c*c+d*d-2*d+1==0)matherror=1;
				a=(1-c*c-d*d)/(c*c+d*d-2*d+1);
				b=(-2*c)/(c*c+d*d-2*d+1);
				if(a==0&&b==0)matherror=1;
			}
			c=log(a*a+b*b)/2;
			d=planeangle(a,b);
			if(d>PI)d-=2*PI;
			if(z==22||z==23)
			{
				a=d;
				b=-c;
			}
			else if(z==24)
			{
				a=-d/2;
				b=c/2;
			}
			else if(z==83)
			{
				a=PI/2+d/2;
				b=-c/2;
			}
			if(numsum==1&&((z==24||z==83)||(fabs(number[1])<=1)))
			{
				optinum=a;
				if(degreemode)optinum*=(180/PI);
			}
			else
			{
				createcomplex(x,a,b);
				if(degreemode)
				{
					data[x].mat[1]*=(180/PI);
					data[x].mat[2]*=(180/PI);
				}
				return;
			}
		}
	}
	else if(z>=57&&z<=59)
	{
		if(numsum>2)matherror=1;
		else if(numsum==2)
		{
			c=number[1];
			d=number[2];
		}
		else
		{
			c=number[1];
			d=0;
		}
		if(numsum<=2)
		{
			if(z==57)
			{
				r=sqrt(sqrt((1+c*c-d*d)*(1+c*c-d*d)+4*c*c*d*d));
				arg=planeangle(1+c*c-d*d,2*c*d)/2;
				if(arg>PI/2)arg-=PI;
				a=r*cos(arg)+c;
				b=r*sin(arg)+d;
			}
			else if(z==58)
			{
				r=sqrt(sqrt((c*c-d*d-1)*(c*c-d*d-1)+4*c*c*d*d));
				arg=planeangle(c*c-d*d-1,2*c*d)/2;
				if(arg>PI/2)arg-=PI;
				a=r*cos(arg)+c;
				b=r*sin(arg)+d;
			}
			else
			{
				if(c*c+d*d-2*c+1==0)matherror=1;
				a=(1-c*c-d*d)/(c*c+d*d-2*c+1);
				b=(2*d)/(c*c+d*d-2*c+1);
				if(a==0&&b==0)matherror=1;
			}
			c=log(a*a+b*b)/2;
			d=planeangle(a,b);
			if(d>PI)d-=2*PI;
			if(z==59)
			{
				c/=2;
				d/=2;
			}
			if(numsum==1&&((z==57)||(z==58&&number[1]>=1)||(z==59&&number[1]>-1&&number[1]<1)))optinum=c;
			else
			{
				createcomplex(x,c,d);
				return;
			}
		}
	}
	else if(z==116)
	{
		double sumabs=0;
		if(sp==0)
		{
			for(i=1;i<=numsum;i++)
				sumabs+=fabs(number[i]);
		}
		else
		{
			if(numsum%2==1)matherror=1;
			for(i=1;i<=numsum;i+=2)
			{
				if(number[i+1]<0)matherror=1;
				sumabs+=fabs(number[i])*number[i+1];
			}
		}
		optinum=sumabs;
	}
	else if(z==118||z==119)
	{
		if(numsum>2)matherror=1;
		else if(numsum==2)
		{
			createcomplex(x,count(number[1],NULL,z,&matherror,degreemode),count(number[2],NULL,z,&matherror,degreemode));
			return;
		}
		else optinum=count(number[1],NULL,z,&matherror,degreemode);
	}
	else if(z==120)
	{
		if(numsum>2)matherror=1;
		else if(numsum==1)
		{
			createcomplex(x,number[1],0);
			return;
		}
		else if(numsum==2)
		{
			if(sp==1)number[2]*=-1;
			createcomplex(x,number[1],number[2]);
			return;
		}
	}
	else if(z==122)
	{
		if(numsum%2==1)matherror=1;
		else
		{
			double sumxy=0,sum2x=0,sumx=0,sum2y=0,sumy=0,avgx=0,avgy=0;
			for(i=1;i<=numsum;i+=2)
			{
				sumxy+=number[i]*number[i+1];
				sumx+=number[i];
				sumy+=number[i+1];
				sum2x+=number[i]*number[i];
				sum2y+=number[i+1]*number[i+1];
			}
			avgx=sumx/(numsum/2);
			avgy=sumy/(numsum/2);
			optinum=(sumxy-(numsum/2)*avgx*avgy)/sqrt((sum2x-(numsum/2)*avgx*avgx)*(sum2y-(numsum/2)*avgy*avgy));
		}
	}
	else if(z==126)
	{
		if((sp<=2&&numsum!=2)||(sp>=3&&numsum!=3))matherror=1;
		else
		{
			if(sp==1)
			{
				data[x].id=x;
				data[x].cla=16;
				data[x].fac=0;
				data[x].r=2;
				data[x].c=1;
				data[x].mat=new double[4];
				data[x].mat[1]=sqrt(number[1]*number[1]+number[2]*number[2]);
				data[x].mat[2]=planeangle(number[1],number[2]);
				if(degreemode)data[x].mat[2]*=(180/PI);
				return;
			}
			else if(sp==2)
			{
				if(degreemode)number[2]/=(180/PI);
				createcomplex(x,number[1]*cos(number[2]),number[1]*sin(number[2]));
				return;
			}
			else
			{
				data[x].id=x;
				data[x].cla=16;
				data[x].fac=0;
				data[x].r=(sp==4||sp==6)?1:3;
				data[x].c=3/data[x].r;
				data[x].mat=new double[5];
				if(sp==3)
				{
					data[x].mat[1]=sqrt(number[1]*number[1]+number[2]*number[2]+number[3]*number[3]);
					data[x].mat[2]=acos(number[3]/data[x].mat[1]);
					data[x].mat[3]=planeangle(number[1],number[2]);
					if(degreemode)
					{
						data[x].mat[2]*=(180/PI);
						data[x].mat[3]*=(180/PI);
					}
				}
				else if(sp==4)
				{
					if(degreemode)
					{
						number[2]/=(180/PI);
						number[3]/=(180/PI);
					}
					data[x].mat[1]=number[1]*sin(number[2])*cos(number[3]);
					data[x].mat[2]=number[1]*sin(number[2])*sin(number[3]);
					data[x].mat[3]=number[1]*cos(number[2]);
				}
				else if(sp==5)
				{
					data[x].mat[1]=sqrt(number[1]*number[1]+number[2]*number[2]);
					data[x].mat[2]=planeangle(number[1],number[2]);
					data[x].mat[3]=number[3];
					if(degreemode)data[x].mat[2]*=(180/PI);
				}
				else if(sp==6)
				{
					if(degreemode)number[2]/=(180/PI);
					data[x].mat[1]=number[1]*cos(number[2]);
					data[x].mat[2]=number[1]*sin(number[2]);
					data[x].mat[3]=number[3];
				}
				else if(sp==7)
				{
					double x0,y0,z0;
					if(degreemode)
					{
						number[2]/=(180/PI);
						number[3]/=(180/PI);
					}
					x0=number[1]*sin(number[2])*cos(number[3]);
					y0=number[1]*sin(number[2])*sin(number[3]);
					z0=number[1]*cos(number[2]);
					data[x].mat[1]=sqrt(x0*x0+y0*y0);
					data[x].mat[2]=planeangle(x0,y0);
					data[x].mat[3]=z0;
					if(degreemode)data[x].mat[2]*=(180/PI);
				}
				else if(sp==8)
				{
					double x0,y0,z0;
					if(degreemode)number[2]/=(180/PI);
					x0=number[1]*cos(number[2]);
					y0=number[1]*sin(number[2]);
					z0=number[3];
					data[x].mat[1]=sqrt(x0*x0+y0*y0+z0*z0);
					data[x].mat[2]=acos(z0/data[x].mat[1]);
					data[x].mat[3]=planeangle(x0,y0);
					if(degreemode)
					{
						data[x].mat[2]*=(180/PI);
						data[x].mat[3]*=(180/PI);
					}
				}
				return;
			}
		}
	}
	else if(z==82)
	{
		if(numsum==1)optinum=count(number[1],sp,82,&matherror,degreemode);
		else
		{
			tempnum=tempnum2=0;
			if(numsum%2==1)matherror=1;
			for(i=1;i<=numsum;i+=2)
			{
				if(number[i+1]<0)matherror=1;
				tempnum+=count(number[i],sp,82,&matherror,degreemode)*number[i+1];
				tempnum2+=number[i+1];
			}
			optinum=round(tempnum/tempnum2,2);
		}
	}
	else if(z==81)
	{
		if(numsum==1)optinum=round(number[1],0);
		else if(numsum==2)optinum=round(number[1],number[2]);
		else if(numsum==3)
		{
			createcomplex(x,round(number[1],number[3]),round(number[2],number[3]));
			return;
		}
		else syntaxerror=1;
	}
	else if(z==127)
	{
		if(numsum%2==1)matherror=1;
		optinum=0;
		for(i=1;i<=numsum;i+=2)
		{
			tempnum=1;
			for(j=1;j<=sp%10;j++)
				tempnum*=number[i];
			for(j=1;j<=sp/10;j++)
				tempnum*=number[i+1];
			optinum+=tempnum;
		}
	}
	else if(z==128)
	{
		bool odd=numsum%2==1;
		int tempsp=sp;
		double pos,*summ;
		if(numsum%2==1)
		{
			pos=number[numsum];
			numsum--;
		}
		if(sp<0)
		{
			if(sp==-1)
			{
				for(i=2;i<=numsum;i+=2)
					number[i]=log(number[i]);
			}
			else if(sp==-2)
			{
				for(i=1;i<=numsum;i+=2)
					number[i]=log(number[i]);
			}
			else if(sp==-3)
			{
				for(i=1;i<=numsum;i++)
					number[i]=log(number[i]);
			}
			else if(sp==-4)
			{
				for(i=1;i<=numsum;i+=2)
					number[i]=1/number[i];
			}
			sp=1;
		}
		summ=new double[2*sp+2];
		data[0].mat=new double[(sp+1)*(sp+1)+2];
		data[0].r=data[0].c=sp+1;
		for(i=0;i<=2*sp;i++)
		{
			tempnum=0;
			for(j=1;j<=numsum;j+=2)
			{
				tempnum2=1;
				for(k=1;k<=i;k++)
					tempnum2*=number[j];
				tempnum+=tempnum2;
			}
			summ[i]=tempnum;
		}
		for(i=1;i<=sp+1;i++)
			for(j=1;j<=sp+1;j++)
				data[0].mat[(i-1)*(sp+1)+j]=summ[i+j-2];
		invmat(0);
		for(i=0;i<=sp;i++)
		{
			tempnum=0;
			for(j=1;j<=numsum;j+=2)
			{
				tempnum2=number[j+1];
				for(k=1;k<=i;k++)
					tempnum2*=number[j];
				tempnum+=tempnum2;
			}
			summ[i]=tempnum;
		}
		data[x].id=x;
		data[x].r=sp+1;
		data[x].c=1;
		data[x].fac=0;
		data[x].cla=16;
		data[x].mat=new double[sp+3];
		j=0;
		for(i=1;i<=sp+1;i++)
		{
			data[x].mat[i]=0;
			for(k=0;k<=sp;k++)
			{
				j++;
				data[x].mat[i]+=data[0].mat[j]*summ[k];
			}
			if(fabs(data[x].mat[i])<1e-11)data[x].mat[i]=0;
		}
		delete []summ;
		if(tempsp<0&&sp==1)
		{
			if(tempsp==-1)
			{
				data[x].mat[1]=exp(data[x].mat[1]);
				data[x].mat[2]=exp(data[x].mat[2]);
			}
			else if(tempsp==-3)data[x].mat[1]=exp(data[x].mat[1]);
		}
		if(odd)
		{
			if(int(pos)>0&&int(pos)<=sp+1)
			{
				optinum=data[x].mat[int(pos)];
				delete []data[x].mat;
			}
			else
			{
				delete []data[x].mat;
				matherr(x);
			}
		}
		else return;
	}
	data[x].id=x;
	data[x].cla=1;
	data[x].num=optinum;
}
/*----------------------------------------------------------------------------------------------------------*/
void createcomplex(const int x,const double a,const double b)
{
	data[x].id=x;
	data[x].cla=16;
	data[x].r=1;
	data[x].c=2;
	data[x].fac=1;
	data[x].mat=new double[4];
	data[x].mat[0]=data[x].mat[3]=0;
	data[x].mat[1]=a;
	data[x].mat[2]=b;
}
/*----------------------------------------------------------------------------------------------------------*/
void calcmatrix(const int x,const int y)
{
	int i,j,k,commasum=0,lastpos=x,numsum=0,semico=0,semico2=0,semico2sum=0,r=1,c=1,rc=1;
	double tempnum,number[2001]={0};
	bool find,csm=1,iscmplx=1;
	data[y].id=y;
	if(data[x].ope==34)iscmplx=0;
	if((data[x].ope>=37&&data[x].ope<=39)||data[x].ope==96||data[x].ope==125)data[x].ope=34;
	for(i=x+1;i<=y;i++)
		if(data[i].id>0&&(data[i].cla==7||data[i].cla==6))
		{
			data[i].id=0;
			calculate(lastpos,i,20);
			find=0;
			for(j=lastpos+1;j<i;j++)
			{
				if(data[j].id>0&&data[j].cla==1)
				{
					data[j].id=0;
					tempnum=data[j].num;
					numsum++;
					if(numsum>=2000)
					{
						matherror=1;
						return;
					}
					number[numsum]=tempnum;
					find=1;
					iscmplx=0;
					break;
				}
				else if(data[j].id>0&&data[j].cla==16)
				{
					data[j].id=0;
					for(k=1;k<=data[j].c*data[j].r;k++)
					{
						numsum++;
						if(numsum>=2000)
						{
							matherror=1;
							delete []data[j].mat;
							return;
						}
						number[numsum]=data[j].mat[k];
						if(data[j].fac!=1)iscmplx=0;
						find=1;
					}
					csm=0;
					if(data[x].ope!=11)
					{
					    if(semico==0)semico=data[j].c;
					    else if(semico!=data[j].c)matherror=1;
					}
					delete []data[j].mat;
					break;
				}
			}
			if(find==0)syntaxerror=1;
			if(data[i].cla==7)commasum++;
			if(data[i].ope==99&&csm)
			{
				csm=0;
				semico=numsum;
			}
			if(data[i].ope==100)
			{
				semico2sum++;
				semico2=numsum;
			}
			lastpos=i;
		}
	dataleft-=(2*commasum);
	if(semico2sum>1||(semico2!=0&&data[x].ope==104))syntaxerror=1;
	else if(semico2!=0)rc=semico2;
	else rc=numsum;
	if(csm||data[x].ope==11)
	{
		r=1;
		c=rc;
	}
	else
	{
		c=semico;
		if(rc%c!=0)matherror=1;
		else r=rc/c;
	}
	if(matherror||syntaxerror)
	{
		matherr(x);
		return;
	}
	if((data[x].ope==11&&rc==1)||(data[x].ope==11&&iscmplx&&rc==2&&fabs(number[2])<1e-15))
	{
		data[x].id=x;
		data[x].cla=1;
		data[x].r=data[x].c=1;
		data[x].num=number[1];
	}
	else
	{
		data[x].id=x;
		data[x].cla=16;
		data[x].fac=iscmplx;
		data[x].r=r;
		data[x].c=c;
		data[x].mat=new double[rc+2];
		data[x].mat[0]=0;
		for(i=1;i<=rc;i++)
			data[x].mat[i]=number[i];
		if(data[x].ope==104)
		{
			data[x].num=rank(r,c,data[x].mat);
			data[x].cla=1;
			delete []data[x].mat;
		}
		else if(data[x].ope==110)
		{
			delete []data[x].mat;
			data[x].mat=new double[4];
			data[x].r=2;
			data[x].c=1;
			data[x].mat[1]=r;
			data[x].mat[2]=c;
		}
		else if(data[x].ope==117)
		{
			if(c-r!=1)
			{
				matherr(x);
				return;
			}
			data[0].r=data[0].c=r;
			data[0].mat=new double[r*r+2];
			j=0;
			for(i=1;i<=r*c;i++)
			{
				if(i%c==0)continue;
				else
				{
					j++;
					data[0].mat[j]=data[x].mat[i];
				}
			}
			invmat(0);
			k=0;
			for(i=1;i<=r;i++)
			{
				data[x].mat[i]=0;
				for(j=c;j<=r*c;j+=c)
				{
					k++;
					data[x].mat[i]+=data[0].mat[k]*data[x].mat[j];
				}
			}
			delete []data[0].mat;
			data[x].r=1;
			data[x].c=r;
			transpose(x);
		}
		else if(data[x].ope==123)
		{
			for(i=1;i<=data[x].r*data[x].c;i++)
				for(j=i+1;j<=data[x].r*data[x].c;j++)
					if(data[x].mat[i]>data[x].mat[j])
					{
						tempnum=data[x].mat[i];
						data[x].mat[i]=data[x].mat[j];
						data[x].mat[j]=tempnum;
					}
		}
		else if(data[x].ope==124)
		{
			for(i=1;i<=data[x].r*data[x].c;i++)
				for(j=i+1;j<=data[x].r*data[x].c;j++)
					if(data[x].mat[i]<data[x].mat[j])
					{
						tempnum=data[x].mat[i];
						data[x].mat[i]=data[x].mat[j];
						data[x].mat[j]=tempnum;
					}
		}
		else if(data[x].ope==72)
		{
			for(i=1;i<=data[x].r*data[x].c;i++)
				data[x].mat[i]=count(data[x].mat[i],NULL,72,&matherror,degreemode);
		}
	}
	if(semico2!=0)
	{
		if(numsum-semico2!=1&&numsum-semico2!=2)syntaxerror=1;
		if(numsum-semico2==1&&!isint(number[semico2+1]))matherror=1;
		if(numsum-semico2==2&&(!isint(number[semico2+1])||!isint(number[semico2+1])))matherror=1;
		if(data[x].cla==1&&numsum-semico2==1&&int(number[semico2+1])!=1)matherror=1;
		if(data[x].cla==1&&numsum-semico2==2&&(int(number[semico2+1])!=1||int(number[semico2+2])!=1))matherror=1;
		if(data[x].cla==16&&numsum-semico2==1&&(int(number[semico2+1])<=0||int(number[semico2+1])>data[x].r*data[x].c))matherror=1;
		if(data[x].cla==16&&numsum-semico2==2&&(int(number[semico2+1])<=0||int(number[semico2+2])<=0||int(number[semico2+1])>data[x].r||int(number[semico2+2])>data[x].c))matherror=1;
		if(syntaxerror||matherror||data[x].cla==1)return;
		if(numsum-semico2==1)i=int(number[semico2+1]);
		else i=(int(number[semico2+1])-1)*data[x].c+int(number[semico2+2]);
		data[x].cla=1;
		data[x].num=data[x].mat[i];
		delete []data[x].mat;
	}
}
/*----------------------------------------------------------------------------------------------------------*/
void calcdiag(const int x,const int y)
{
	double tempnum,number[4001]={0},**a;
	int i,j,k,l,m,semicosum=0,semico=0,commasum=0,semico2=0,semico2sum=0,lastpos=x,numsum=0,numeach[200]={0},r=0,tempr=0,tempc=0,tempr2=0,tempc2=0;
	bool find;
	data[y].id=y;
	for(i=x+1;i<=y;i++)
		if(data[i].id>0&&(data[i].cla==7||data[i].cla==6))
		{
			data[i].id=0;
			calculate(lastpos,i,20);
			find=0;
			for(j=lastpos+1;j<i;j++)
			{
				if(data[j].id>0&&data[j].cla==1)
				{
					data[j].id=0;
					tempnum=data[j].num;
					numsum++;
					if(numsum>=4000)
					{
						matherror=1;
						return;
					}
					number[numsum]=tempnum;
					find=1;
					break;
				}
				else if(data[j].id>0&&data[j].cla==16)
				{
					data[j].id=0;
					for(k=1;k<=data[j].c*data[j].r;k++)
					{
						numsum++;
						if(numsum>=4000)
						{
							matherror=1;
							delete []data[j].mat;
							return;
						}
						number[numsum]=data[j].mat[k];
						find=1;
					}
					delete []data[j].mat;
					break;
				}
			}
			if(find==0)syntaxerror=1;
			if(data[i].cla==7)commasum++;
			if(data[i].ope==99||(data[i].cla==6&&semico!=0))
			{
				semicosum++;
				numeach[semicosum]=numsum-semico;
				semico=numsum;
			}
			if(data[i].ope==100)
			{
				semico2sum++;
				semico2=numsum;
			}
			lastpos=i;
		}
	dataleft-=(2*commasum);
	if(data[x].ope==105)
	{
		if(semico2!=0)
		{
			matherr(x);
			return;
		}
		if(semico==0)
		{
			data[x].id=x;
			data[x].cla=16;
			data[x].fac=0;
			data[x].mat=new double[numsum*numsum+2];
			data[x].r=data[x].c=numsum;
			for(i=0;i<=numsum*numsum;i++)
				data[x].mat[i]=0;
			j=0;
			for(i=1;i<=numsum*numsum;i+=numsum+1)
			{
				j++;
				data[x].mat[i]=number[j];
			}
		}
		else
		{
			for(i=1;i<=semicosum;i++)
			{
				r+=int(gauss(sqrt(numeach[i])));
				if(!isint(sqrt(numeach[i])))
				{
					matherr(x);
					return;
				}
			}
			data[x].id=x;
			data[x].cla=16;
			data[x].fac=0;
			data[x].r=data[x].c=r;
			data[x].mat=new double[r*r+2];
			a=new double*[r+2];
			for(i=0;i<=r+1;i++)
		        a[i]=new double[r+2];
			k=m=0;
			for(i=1;i<=r;i++)
				for(j=1;j<=r;j++)
				{
					k++;
				    data[x].mat[k]=a[i][j]=0;
				}
			tempr=tempc=0;
			for(i=1;i<=semicosum;i++)
			{
				l=int(gauss(sqrt(numeach[i])));
				for(j=1;j<=l;j++)
					for(k=1;k<=l;k++)
					{
						m++;
						a[tempr+j][tempc+k]=number[m];
					}
				tempr+=l;
				tempc+=l;
			}
			k=0;
			for(i=1;i<=r;i++)
				for(j=1;j<=r;j++)
				{
					k++;
					data[x].mat[k]=a[i][j];
				}
				for(i=0;i<r+1;i++)
		            delete []a[i];
	        delete []a;
		}
	}
	else if(data[x].ope==98)
	{
		if((semico2!=0&&numsum-semico2!=2)||semico2sum>1)
		{
			syntaxerror=1;
			matherr(x);
			return;
		}
		else if(semico2!=0)
		{
			numsum-=2;
			if(semico!=0)
			{
				numeach[semicosum]-=2;
			}
		}
		if(semico==0)
		{
			if(numsum==1)
			{
				data[x].id=x;
				data[x].cla=16;
				data[x].fac=0;
				data[x].r=data[x].c=1;
				data[x].mat=new double[3];
				data[x].mat[1]=number[1];
			}
			else
			{
				if(number[1]*number[2]>2000||number[1]<=0||number[2]<=0||!isint(number[1])||!isint(number[2]))
				{
					matherr(x);
					return;
				}
				data[x].id=x;
				data[x].cla=16;
				data[x].fac=0;
				data[x].r=number[1];
				data[x].c=number[2];
				data[x].mat=new double[data[x].r*data[x].c+2];
				if(numsum==2)
				{
					for(i=1;i<=data[x].r*data[x].c;i++)
						data[x].mat[i]=0;
				}
				else if(numsum==3)
				{
					for(i=1;i<=data[x].r*data[x].c;i++)
						data[x].mat[i]=number[3];
				}
				else
				{
					j=2;
					for(i=1;i<=data[x].r*data[x].c;i++)
					{
						j++;
						data[x].mat[i]=j>numsum?0:number[j];
					}
				}
			}
		}
		else
		{
			if(numeach[1]==1||number[1]*number[2]>2000||number[1]<=0||number[2]<=0||!isint(number[1])||!isint(number[2]))
			{
				matherr(x);
				return;
			}
			for(i=2;i<=semicosum;i++)
				if(numeach[i]==1||numeach[i]==2)
				{
					matherr(x);
					return;
				}
			data[x].id=x;
			data[x].cla=16;
			data[x].fac=0;
			data[x].r=number[1];
			data[x].c=number[2];
			data[x].mat=new double[data[x].r*data[x].c+2];
			a=new double*[data[x].r+2];
			for(i=0;i<=data[x].r+1;i++)
		        a[i]=new double[data[x].c+2];
			k=2;
			for(i=1;i<=data[x].r;i++)
				for(j=1;j<=data[x].c;j++)
				{
					if(numeach[1]==3)a[i][j]=number[3];
					else
					{
						k++;
					    a[i][j]=k>numeach[1]?0:number[k];
					}
				}
			k=numeach[1];
			for(i=2;i<=semicosum;i++)
			{
				if(numeach[i]==3)
				{
					k+=3;
					if(!isint(number[k-2])||!isint(number[k-1])||number[k-2]>data[x].r||number[k-1]>data[x].c||number[k-2]<=0||number[k-1]<=0)
					{
						matherr(x);
						return;
					}
					tempr=number[k-2];
					tempc=number[k-1];
					a[tempr][tempc]=number[k];
				}
				else
				{
					lastpos=k;
					k+=4;
					if(!isint(number[k-3])||!isint(number[k-2])||!isint(number[k-1])||!isint(number[k])||number[k-3]>data[x].r||number[k-3]<=0||number[k-2]>data[x].c||number[k-2]<=0||number[k-1]>data[x].r||number[k-1]<=0||number[k]>data[x].c||number[k]<=0)
					{
						matherr(x);
						return;
					}
					tempr=number[k-3];
					tempc=number[k-2];
					tempr2=number[k-1];
					tempc2=number[k];
					j=4;
					for(l=tempr;l<=tempr2;l++)
						for(m=tempc;m<=tempc2;m++)
						{
							if(numeach[i]==5)a[l][m]=number[k+1];
							else if(j<numeach[i])
							{
								j++;
								k++;
								a[l][m]=number[k];
							}
							else break;
						}
					k=lastpos+numeach[i];
				}
			}
			k=0;
			for(i=1;i<=data[x].r;i++)
				for(j=1;j<=data[x].c;j++)
				{
					k++;
					data[x].mat[k]=a[i][j];
				}
				for(i=0;i<r+1;i++)
		            delete []a[i];
	        delete []a;
		}
		if(semico2!=0)
		{
			if(!isint(number[numsum+1])||!isint(number[numsum+2])||number[numsum+1]>data[x].r||number[numsum+1]<=0||number[numsum+2]>data[x].c||number[numsum+2]<=0)
			{
				matherror=1;
				number[numsum+1]=number[numsum+2]=0;
			}
			tempr=number[numsum+1];
			tempc=number[numsum+2];
			i=(tempr-1)*data[x].c+tempc;
			data[x].num=data[x].mat[i];
			data[x].cla=1;
			delete []data[x].mat;
		}
	}
}
/*----------------------------------------------------------------------------------------------------------*/
inline void matherr(const int x)
{
	matherror=1;
	data[x].id=x;
	data[x].cla=1;
	data[x].num=0;
}
/*----------------------------------------------------------------------------------------------------------*/
void detcm(const int x,const int n,const double number[ ])
{
	double **re=new double*[n+2],**im=new double*[n+2],bre,bim,temp,tempre,tempim;
	int i,j=0,k=1,m=1;
	bool t;
	for(i=0;i<=n+1;i++)
	{
		re[i]=new double[n+2];
		im[i]=new double[n+2];
	}
	for(i=1;i<=2*n*n;i+=2)
	{
		j++;
		if(j>n)
		{
			k++;
			j=1;
		}
		re[k][j]=number[i];
		im[k][j]=number[i+1];
	}
	for(i=1;i<=n;i++)
	{
		t=(re[i][i]!=0||im[i][i]!=0);
		if(!t)
		{
			for(j=i+1;j<=n;j++)
				if(re[j][i]!=0||im[j][i]!=0)
				{
					t=1;
					for(k=1;k<=n;k++)
					{
						temp=re[i][k];
						re[i][k]=re[j][k];
						re[j][k]=temp;
						temp=im[i][k];
						im[i][k]=im[j][k];
						im[j][k]=temp;
					}
					m*=-1;
					break;
				}
			if(!t)
			{
				for(i=0;i<=n+1;i++)
				{
					delete []re[i];
					delete []im[i];
				}
				delete []re;
				delete []im;
				createcomplex(x,0,0);
				return;
			}
		}
		for(j=i+1;j<=n;j++)
		{
			bre=-(re[j][i]*re[i][i]+im[i][i]*im[j][i])/(re[i][i]*re[i][i]+im[i][i]*im[i][i]);
			bim=(re[j][i]*im[i][i]-im[j][i]*re[i][i])/(re[i][i]*re[i][i]+im[i][i]*im[i][i]);
			for(k=i;k<=n;k++)
			{
				tempre=bre*re[i][k]-bim*im[i][k];
				tempim=bim*re[i][k]+bre*im[i][k];
				re[j][k]+=tempre;
				im[j][k]+=tempim;
			}
		}
	}
	bre=m;
	bim=0;
	for(i=1;i<=n;i++)
	{
		temp=bre*re[i][i]-bim*im[i][i];
		bim=bre*im[i][i]+bim*re[i][i];
		bre=temp;
	}
	for(i=0;i<=n+1;i++)
	{
		delete []re[i];
		delete []im[i];
	}
	delete []re;
	delete []im;
	createcomplex(x,bre,bim);
}
/*----------------------------------------------------------------------------------------------------------*/
void transpose(const int x)
{
	if(data[x].cla==1)return;
	int r,c,i,j,k=0;
	double a[2001]={0};
	r=data[x].r;
	c=data[x].c;
	for(i=1;i<=c;i++)
		for(j=1;j<=r;j++)
		{
			k++;
			a[k]=data[x].mat[(j-1)*c+i];
		}
	for(i=1;i<=r*c;i++)
		data[x].mat[i]=a[i];
	data[x].r=c;
	data[x].c=r;
	data[x].fac=0;
}
/*----------------------------------------------------------------------------------------------------------*/
void adjmat(const int x)
{
	if(data[x].cla==1)
	{
		data[x].num=1;
		return;
	}
	int i,j,k,m,r,c,r1,c1;
	double a[2001],b[2001];
	if(data[x].r!=data[x].c)
	{
		matherror=1;
		return;
	}
	for(i=1;i<=data[x].r*data[x].r;i++)
	{
	    k=0;
	    r=i%data[x].r!=0?i/data[x].r+1:i/data[x].r;
	    c=i%data[x].r!=0?i%data[x].r:data[x].r;
	    m=abs(r-c)%2==0?1:-1;
	    for(j=1;j<=data[x].r*data[x].r;j++)
		{
		    r1=j%data[x].r!=0?j/data[x].r+1:j/data[x].r;
	        c1=j%data[x].r!=0?j%data[x].r:data[x].r;
		    if(r!=r1&&c!=c1)
			{
			    k++;
			    b[k]=data[x].mat[j];
			}
		}
	    a[i]=m*det(data[x].r-1,b);
	}
    for(i=1;i<=data[x].r*data[x].r;i++)
	    data[x].mat[i]=a[i];
	transpose(x);
}
/*----------------------------------------------------------------------------------------------------------*/
void invmat(const int x)
{
	int i;
	double temp;
	temp=det(data[x].r,data[x].mat);
	if(temp==0)matherror=1;
	adjmat(x);
	for(i=1;i<=data[x].r*data[x].r;i++)
		data[x].mat[i]/=temp;
}
/*----------------------------------------------------------------------------------------------------------*/
inline void oneomatrix(const int x)
{
	data[x].cla=1;
	data[x].num=data[x].mat[1];
	delete []data[x].mat;
}
/*----------------------------------------------------------------------------------------------------------*/
void calcordinfun(const int x,const int y)
{
	int i,j,k,varsum=0,prev=x,next=x,z=data[x].ope-64,dts=datasum,dtl,length,ansr,ansc,ansf;
	bool fmd=funmode,ansism=0,ansisd=0;
	double funvari[6]={0},funvarmem[5]={0},ansv,ansm[2001]={0};
	char in[5000];
	countfeverexpressiondata tempcfr[3000];
	iterative++;
	if(iterative>=50)
	{
		matherror=1;
		for(i=x+1;i<y;i++)
			data[i].id=0;
		dataleft=1;
		data[x].id=x;
		data[x].cla=1;
		data[x].num=0;
		iterative--;
		return;
	}
	for(i=0;i<5000;i++)
		in[i]=input[i];
	for(i=1;i<=4;i++)
		funvarmem[i]=var[i].temp;
	data[y].id=y;
	length=strlen(ordinfun[z]);
	dataleft++;
	for(i=x+1;i<=y;i++)
	{
		if(data[i].id>0&&(data[i].cla==7||data[i].cla==6))
		{
			prev=next;
			next=i;
			data[i].id=0;
			dataleft--;
			calculate(prev,next,20);
			for(j=prev+1;j<next;j++)
			{
				if(data[j].id>0&&data[j].cla==1)
				{
					varsum++;
					if(varsum>5)
					{
						varsum=5;
						matherror=1;
					}
					funvari[varsum]=data[j].num;
					data[j].id=0;
					dataleft--;
					break;
				}
				else if(data[j].id>0&&data[j].cla==16)
				{
					for(k=1;k<=data[j].r*data[j].c;k++)
					{
						varsum++;
						if(varsum>5)
						{
							varsum=5;
							matherror=1;
						}
						funvari[varsum]=data[j].mat[k];
					}
					dataleft--;
					data[j].id=0;
					delete []data[j].mat;
					break;
				}
			}
		}
	}
	dtl=dataleft;
	for(i=1;i<=dts;i++)
		tempcfr[i]=data[i];
	j=0;
	for(i=0; ;i++)
	{
		if(ordinfun[z][i]=='=')break;
		if(ordinfun[z][i]>=87&&ordinfun[z][i]<=90)
		{
			j++;
			k=ordinfun[z][i]-87>0?ordinfun[z][i]-87:4;
			var[k].temp=funvari[j];
		}
	}
	if(varsum!=j)matherror=1;
	for(j=0;j<5000;j++)
		input[j]='\0';
	for(j=i+1;j<length;j++)
		input[j-i-1]=ordinfun[z][j];
	translate();
	dataleft=datasum;
	maincalc();
	for(i=1;i<=datasum;i++)
		if(data[i].id>0)
		{
			if(data[i].cla==16)
			{
				ansism=1;
				ansr=data[i].r;
				ansc=data[i].c;
				ansf=data[i].fac;
				for(j=1;j<=ansr*ansc;j++)
					ansm[j]=data[i].mat[j];
				delete []data[i].mat;
			}
			else if(data[i].cla==17)
			{
				ansisd=1;
				ansf=data[i].fac;
				ansr=data[i].r;
				ansc=data[i].c;
			}
			ansv=data[i].num;
			break;
		}
	datasum=dts;
	dataleft=dtl+1;
	funmode=fmd;
	for(i=0;i<5000;i++)
		input[i]=in[i];
	for(i=1;i<=4;i++)
		var[i].temp=funvarmem[i];
	for(i=1;i<=dts;i++)
		data[i]=tempcfr[i];
	data[x].id=x;
	data[x].cla=1;
	if(!ansism&&!ansisd)data[x].num=ansv;
	else if(ansisd)
	{
		data[x].cla=17;
		data[x].fac=ansf;
		data[x].r=ansr;
		data[x].c=ansc;
	}
	else
	{
		data[x].cla=16;
		data[x].fac=ansf;
		data[x].r=ansr;
		data[x].c=ansc;
		data[x].mat=new double[ansr*ansc+2];
		for(i=1;i<=ansr*ansc;i++)
			data[x].mat[i]=ansm[i];
	}
	iterative--;
}
/*----------------------------------------------------------------------------------------------------------*/
int maincalc()
{
	int i,k,wtime=0,temppos=0;
	if(dataleft==1&&data[1].cla>=9&&data[1].cla<=12)
	{
		data[1].num=var[data[1].cla-8].temp;
		data[1].cla=1;
		return 0;
	}
	while(dataleft>1&&wtime<400)
	{
        wtime++;
		leftbrac=rightbrac=0;
        k=searchbrac();
   	    if(leftbrac==0&&rightbrac==0)
		{
			rightbrac=datasum+1;
		}
		if(k==43||k==44||k==73||k==74||k==31||k==32||k==12||k==15||k==19||(k>=84&&k<=94)||k==97||k==106||(k>=114&&k<=116)||(ismatrix==1&&k!=11&&k!=34&&k!=125&&k!=72)||(k>=22&&k<=24)||k==83||(k>=57&&k<=59)||(k>=120&&k<=122)||(k>=126&&k<=128)||k==82||k==81)calcvector(leftbrac,rightbrac);
		else if(k==104||k==110||k==117||k==123||k==124||ismatrix)calcmatrix(leftbrac,rightbrac);
		else if(k==98||k==105)calcdiag(leftbrac,rightbrac);
		else if(k>=65&&k<=70)calcordinfun(leftbrac,rightbrac);
		else
		{
			calculate(leftbrac,rightbrac,20);
			temppos=0;
			for(i=leftbrac+1;i<rightbrac;i++)
			{
				if(data[i].id>0&&(data[i].cla!=16||(data[i].cla==16&&data[i].r==1&&data[i].c==1)||leftbrac==0))
				{
					if(data[i].cla==16&&data[i].r==1&&data[i].c==1)oneomatrix(i);
					data[i].num=count(data[i].num,0,k,&matherror,degreemode);
					temppos=i;
			        break;
				}
				else if(data[i].id>0&&data[i].cla==16)
				{
					temppos=i;
					break;
				}
			}
			if(temppos==0)matherror=1;
			if(rightbrac<=datasum&&data[rightbrac+1].cla==15)
			{
				if(fabs(data[temppos].num)<1e-15)
				{
				    ondisplay=0;
				    dataleft=1;
				    for(i=rightbrac;i<=datasum;i++)
					    data[i].id=0;
				    return wtime;
				}
				else
				{
					dataleft-=2;
					data[temppos].id=data[rightbrac+1].id=0;
				}
			}
		}
	}
	for(i=1;i<=datasum;i++)
		if(data[i].id>0)
		{
			if(data[i].cla>=9&&data[i].cla<=12)
			{
				data[i].num=var[data[i].cla-8].temp;
				data[i].cla=1;
			}
			break;
		}
	return wtime;
}
/*----------------------------------------------------------------------------------------------------------*/
void calculate(const int x,const int y,const int z)
{
	int i,j,t,k,opesum,opelist[10]={0},m1p,m2p;
	bool find,mat1=0,mat2=0,date1=0,date2=0;
	double temp1,temp2;
	if(z==20)
	{
		tranfunvar(x,y);
		for(i=x+1;i<y;i++)
			if(data[i].id>0&&data[i].cla==3)
			{
				data[i].id=0;
				dataleft--;
				find=0;
				for(j=i-1;j>x;j--)
					if(data[j].id>0&&(data[j].cla==1||data[j].cla==16||data[j].cla==17))
					{
						find=1;
						break;
					}
				if(find==0||data[j].cla==17)syntaxerror=1;
				if(data[j].cla==1||data[j].cla==16)
				{
				    if(data[i].ope==8&&data[i+1].id>0&&data[i+1].cla==3&&data[i+1].ope==8)
					{
					    data[i+1].id=0;
					    dataleft--;
						if(data[j].r==1&&data[j].c==2&&data[j].cla==16&&data[j].fac==1&&fabs(data[j].mat[2])<1e-15)
						{
							data[j].num=data[j].mat[1];
							data[j].cla=1;
							delete []data[j].mat;
						}
					    data[j].num=factdouble(data[j].num,&matherror);
					}
				    else if(data[i].ope==8)
					{
						if(data[j].r==1&&data[j].c==2&&data[j].cla==16&&data[j].fac==1&&fabs(data[j].mat[2])<1e-15)
						{
							data[j].num=data[j].mat[1];
							data[j].cla=1;
							delete []data[j].mat;
						}
						data[j].num=fact(data[j].num,&matherror);
					}
				    else if(data[i].ope==40)data[j].num/=100;
				    else if(data[i].ope==41)
					{
						if(data[j].cla==16)matrixpower(j,2);
						else data[j].num*=data[j].num;
					}
				    else if(data[i].ope==42)
					{
						if(data[j].cla==1)data[j].num*=(data[j].num*data[j].num);
						else matrixpower(j,3);
					}
					else if(data[i].ope==101)transpose(j);
					else if(data[i].ope==102)adjmat(j);
					else matherror=1;
				}
			}
		for(i=x+1;i<y;i++)
			if(data[i].id>0&&data[i].cla==13)
			{
				data[i].id=0;
				dataleft--;
				find=0;
				for(j=i+1;j<y;j++)
					if(data[j].id>0&&(data[j].cla==1||data[j].cla==16||data[j].cla==17))
					{
						find=1;
						break;
					}
				if(find==0||data[j].cla==17)syntaxerror=1;
				if(data[i].ope==48&&data[j].cla==1)data[j].num=loginot(data[j].num);
				else if(data[i].ope==48&&data[j].cla==16)
				{
					for(k=1;k<=data[j].r*data[j].c;k++)
						data[j].mat[k]=loginot(data[j].mat[k]);
				}
				else if(data[i].ope==78&&data[j].cla==1)data[j].num=~(int(data[j].num));
				else if(data[i].ope==78&&data[j].cla==16)
				{
					for(k=1;k<=data[j].r*data[j].c;k++)
						data[j].mat[k]=~(int(data[j].mat[k]));
				}
				else matherror=1;
			}
		for(i=x+1;i<y;i++)
			if(data[i].id>0&&data[i].cla==4)
			{
				data[i].id=0;
				dataleft--;
				find=0;
				if(data[i].ope==10)
				{
					for(j=i+1;j<y;j++)
						if(data[j].id>0&&(data[j].cla==1||data[j].cla==4||data[j].cla==16||data[j].cla==17))
						{
							find=1;
							if(data[j].cla==4&&data[j].ope==10)data[j].ope=9;
							else if(data[j].cla==4&&data[j].ope==9)data[j].ope=10;
							else if(data[j].cla==1)data[j].num*=-1;
							else if(data[j].cla==16)
							{
								for(k=1;k<=data[j].r*data[j].c;k++)
									data[j].mat[k]*=-1;
							}
							else if(data[j].cla==17)data[j].fac*=-1;
							break;
						}
					if(find==0)syntaxerror=1;
				}
			}
		calculate(x,y,9);
	}
	else
	{
		t=0;
		for(i=x+1;i<y;i++)
			if(data[i].id>0&&(data[i].cla==2||data[i].cla==14||data[i].cla==15))
			{
				t=1;
				break;
			}
		if(t==0)return;
		if(z==9)
		{
			opesum=5;
			opelist[1]=107;
			opelist[2]=108;
			opelist[3]=109;
			opelist[4]=111;
			opelist[5]=112;
		}
		else if(z==8)
		{
			opesum=1;
			opelist[1]=121;
		}
		else if(z==7)
		{
			opesum=2;
			opelist[1]=29;
			opelist[2]=30;
		}
		else if(z==6)
		{
			opesum=2;
			opelist[1]=79;
			opelist[2]=80;
		}
		else if(z==5)
		{
			opesum=3;
			opelist[1]=75;
			opelist[2]=76;
			opelist[3]=77;
		}
		else if(z==4)
		{
			opesum=4;
			opelist[1]=45;
			opelist[2]=46;
			opelist[3]=47;
			opelist[4]=95;
		}
		else if(z==3)
		{
			opesum=1;
			opelist[1]=5;
		}
		else if(z==2)
		{
			opesum=4;
			opelist[1]=3;
			opelist[2]=4;
			opelist[3]=6;
			opelist[4]=7;
		}
		else if(z==1)
		{
			opesum=2;
			opelist[1]=1;
			opelist[2]=2;
		}
		for(i=x+1;i<y;i++)
		{
			mat1=mat2=date1=date2=0;
			if(data[i].id>0&&data[i].cla==2&&include(data[i].ope,z))
			{
				find=0;
				for(j=i+1;j<y;j++)
					if(data[j].id>0&&(data[j].cla==1||data[j].cla==16||data[j].cla==17))
					{
						if(data[j].cla==1)temp2=data[j].num;
						else if(data[j].cla==16)mat2=1;
						else if(z==1)date2=1;
						m2p=j;
						data[j].id=0;
						find=1;
						break;
					}
				if(find==0||(z!=1&&data[j].cla==17))syntaxerror=1;
				find=0;
				for(j=i-1;j>x;j--)
					if(data[j].id>0&&(data[j].cla==1||data[j].cla==16||data[j].cla==17))
					{
						if(data[j].cla==1)temp1=data[j].num;
						else if(data[j].cla==16)mat1=1;
						else if(z==1)date1=1;
						m1p=j;
						find=1;
						break;
					}
				if(find==0||(z!=1&&data[j].cla==17))syntaxerror=1;
				dataleft-=2;
				if(z==8||z==9)countmatrix(m1p,m2p,data[i].ope);
				else if(z==1&&date1&&date2)
				{
					data[j].cla=1;
					if(data[i].ope==1)data[m2p].fac*=-1;
					data[j].num=daysbetween(data[m1p].fac,data[m1p].r,data[m1p].c,data[m2p].fac,data[m2p].r,data[m2p].c,&matherror);
				}
				else if(z==1&&date1)
				{
					if(mat2)
					{
						if((data[m2p].fac==1&&data[m2p].r==1&&data[m2p].c==2&&fabs(data[m2p].mat[2])<1e-15)||(data[m2p].r==1&&data[m2p].c==1))oneomatrix(m2p);
						else
						{
							matherr(mat2);
							syntaxerror=1;
						}
					}
					data[m2p].id=0;
					if(data[i].ope==2)data[m2p].num*=-1;
					countdate(&data[m1p].fac,&data[m1p].r,&data[m1p].c,data[m2p].num,&matherror);
				}
				else if(z==1&&date2)
				{
					if(mat1)
					{
						if((data[m1p].fac==1&&data[m1p].r==1&&data[m1p].c==2&&fabs(data[m1p].mat[2])<1e-15)||(data[m1p].r==1&&data[m1p].c==1))oneomatrix(m1p);
						else
						{
							matherr(mat1);
							syntaxerror=1;
						}
					}
					data[j].id=0;
					data[m2p].id=m2p;
					if(data[i].ope==2)data[m2p].fac*=-1;
					countdate(&data[m2p].fac,&data[m2p].r,&data[m2p].c,data[m2p].num,&matherror);
				}
				else if(!mat1&&!mat2)
				{
					if(z!=3)data[j].num=count(temp1,temp2,data[i].ope,&matherror,degreemode);
					else complexpower(m1p,m2p);
				}
				else if(mat1&&!mat2)
				{
					if(data[i].ope==1||data[i].ope==2)
					{
						if(data[m1p].r==1&&data[m1p].c==2)
						{
							data[m1p].fac=1;
							data[m1p].mat[1]=count(data[m1p].mat[1],data[m2p].num,data[i].ope,&matherror,degreemode);
						}
						else if(data[m1p].r==data[m1p].c)
						{
							for(j=1;j<=data[m1p].r*data[m1p].r;j+=data[m1p].r+1)
								data[m1p].mat[j]=count(data[m1p].mat[j],data[m2p].num,data[i].ope,&matherror,degreemode);
						}
						else matherror=1;
					}
					else if(data[i].ope==3||data[i].ope==4||data[i].ope==6||data[i].ope==7||(data[i].ope==103&&data[m1p].r==1&&data[m1p].c==2&&data[m1p].fac==1)||z==6)
					{
						for(j=1;j<=data[m1p].r*data[m1p].c;j++)
							data[m1p].mat[j]=count(data[m1p].mat[j],data[m2p].num,data[i].ope,&matherror,degreemode);
					}
					else if(data[i].ope==5)matrixpower(m1p,data[m2p].num);
					else if(z==5&&data[m1p].r==1&&data[m1p].c==2&&data[m1p].fac==1&&fabs(data[m1p].mat[2])<1e-15)
					{
						data[m1p].num=count(data[m1p].mat[1],data[m2p].num,data[i].ope,&matherror,degreemode);
						data[m1p].cla=1;
						delete []data[m1p].mat;
					}
					else matherror=1;
				}
				else if(mat2&&!mat1)
				{
					if(data[i].ope==3||data[i].ope==1||data[i].ope==2||(data[i].ope==103&&data[m1p].r==1&&data[m1p].c==2&&data[m1p].fac==1))
					{
					    data[m2p].id=m2p;
					    data[m1p].id=0;
						if(data[i].ope==1||data[i].ope==2)
						{
							if(data[m2p].r==1&&data[m2p].c==2)
							{
								data[m2p].fac=1;
								data[m2p].mat[1]=count(data[m1p].num,data[m2p].mat[1],data[i].ope,&matherror,degreemode);
								data[m2p].mat[2]=count(0,data[m2p].mat[2],data[i].ope,&matherror,degreemode);
							}
							else if(data[m2p].r==data[m2p].c)
							{
								if(data[i].ope==2)
								{
									for(j=1;j<=data[m2p].r*data[m2p].c;j++)
										data[m2p].mat[j]*=-1;
								}
								for(j=1;j<=data[m2p].r*data[m2p].r;j+=data[m2p].r+1)
									data[m2p].mat[j]+=data[m1p].num;
							}
							else matherror=1;
						}
						else
						{
							for(j=1;j<=data[m2p].r*data[m2p].c;j++)
								data[m2p].mat[j]*=data[m1p].num;
						}
					}
					else if(data[i].ope==4)
					{
						data[m2p].id=m2p;
					    data[m1p].id=0;
						if(data[m2p].r==1&&data[m2p].c==2&&data[m2p].mat[1]==0&&data[m2p].mat[2]==0)
						{
							delete []data[m2p].mat;
							matherr(m2p);
						}
						else if(data[m2p].r==1&&data[m2p].c==2)
						{
							data[m2p].fac=1;
							temp1=data[m1p].num*data[m2p].mat[1]/(data[m2p].mat[1]*data[m2p].mat[1]+data[m2p].mat[2]*data[m2p].mat[2]);
							temp2=-data[m1p].num*data[m2p].mat[2]/(data[m2p].mat[1]*data[m2p].mat[1]+data[m2p].mat[2]*data[m2p].mat[2]);
							data[m2p].mat[1]=temp1;
							data[m2p].mat[2]=temp2;
						}
						else
						{
						    invmat(m2p);
						    for(j=1;j<=data[m2p].r*data[m2p].c;j++)
							    data[m2p].mat[j]*=data[m1p].num;
						}
					}
					else if(z==5&&data[m2p].r==1&&data[m2p].c==2&&data[m2p].fac==1&&fabs(data[m2p].mat[2])<1e-15)
					{
						data[m2p].id=m2p;
						data[m1p].id=0;
						data[m2p].num=count(data[m2p].mat[1],data[m1p].num,data[i].ope,&matherror,degreemode);
						data[m2p].cla=1;
						delete []data[m2p].mat;
					}
					else if(data[i].ope==5)complexpower(m1p,m2p);
					else matherror=1;
				}
				else countmatrix(m1p,m2p,data[i].ope);
				data[i].id=0;
			}
		}
		if(z>1)calculate(x,y,z-1);
	}
	if(z==1)
	{
		int logicla=0,n1p=0,cmpsum=0,logi[6]={0},str,stc;
		double cmplist[2001]={0};
		mat1=0;
		for(i=x+1;i<y;i++)
		{
			if(data[i].id>0&&data[i].cla==16&&((data[i].r==1&&data[i].c==1)||(data[i].r==1&&data[i].c==2&&data[i].fac==1&&fabs(data[i].mat[2])<1e-15)))oneomatrix(i);
			if(data[i].id>0&&data[i].cla==16)mat1=1;
			if(data[i].id>0&&data[i].cla==14)logi[data[i].ope-49]++;
		}
		if(mat1&&(logi[0]+logi[1]+logi[2]+logi[3]+logi[4]+logi[5])>0)
		{
			if(logi[0]+logi[1]+logi[3]+logi[4]+logi[5]>0)matherror=1;
			mat1=0;
			for(i=x+1;i<y;i++)
			{
				if(data[i].id>0&&data[i].cla==14)
				{
					dataleft--;
					data[i].id=0;
				}
				if(data[i].id>0&&data[i].cla==16)
				{
					if(!mat1)
					{
						str=data[i].r;
						stc=data[i].c;
						mat1=1;
						for(j=1;j<=str*stc;j++)
							cmplist[j]=data[i].mat[j];
						delete []data[i].mat;
						data[i].cla=1;
						data[i].num=1;
						n1p=i;
					}
					else
					{
						dataleft--;
						data[i].id=0;
						if(str!=data[i].r||stc!=data[i].c)data[n1p].num=0;
						for(j=1;j<=data[i].r*data[i].c&&data[n1p].num==1;j++)
							if(fabs(cmplist[j]-data[i].mat[j])>1e-15)
							{
								data[n1p].num=0;
								break;
							}
						data[i].cla=1;
						delete []data[i].mat;
					}
				}
				if(data[i].id>0&&(data[i].cla==1||data[i].cla==17))
				{
					dataleft--;
					data[i].id=0;
					matherror=1;
				}
			}
			return;
		}
		if(((logi[0]+logi[3])*(logi[1]+logi[4]))>0)matherror=1;
		if(((logi[0]+logi[1]+logi[2]+logi[3]+logi[4])*logi[5])>0)matherror=1;
		if((logi[0]+logi[1]+logi[2]+logi[3]+logi[4]+logi[5])>0)
		{
			if(matherror==0&&logi[5]>0)
			{
				t=1;
				for(i=x+1;i<y;i++)
				{
					if(data[i].id>0&&data[i].cla==15&&data[i].ope==60)break;
					if(data[i].id>0&&data[i].cla==14)
					{
						dataleft--;
						data[i].id=0;
					}
					else if(data[i].id>0&&data[i].cla==1)
					{
						if(cmpsum>0)
						{
							dataleft--;
						    data[i].id=0;
						}
						else n1p=i;
						if(t==1)
						{
						    for(j=1;j<=cmpsum;j++)
							    if(fabs(cmplist[j]-data[i].num)<1e-15)
								{
									t=0;
									break;
								}
							if(t==1)
							{
								cmpsum++;
								cmplist[cmpsum]=data[i].num;
							}
						}
					}
				}
				data[n1p].num=t;
			}
			else if(matherror==0)
			{
				t=1;
				for(i=x+1;i<y;i++)
				{
					if(data[i].id>0&&data[i].cla==15&&data[i].ope==60)break;
					if(data[i].id>0&&data[i].cla==1)
					{
						if(cmpsum==0)n1p=i;
						else
						{
							dataleft--;
							data[i].id=0;
						}
						cmpsum++;
						cmplist[cmpsum]=data[i].num;
					}
				}
				j=0;
				for(i=x+1;i<y;i++)
					if(data[i].id>0&&data[i].cla==14)
					{
						j++;
						dataleft--;
						data[i].id=0;
						if(t==1)t=cmp(cmplist[j],cmplist[j+1],data[i].ope);
					}
				data[n1p].num=t;
			}
			else
			{
				for(i=x+1;i<y;i++)
					if(data[i].id>0&&(data[i].cla==1||data[i].cla==14))
					{
						if(data[i].cla==1&&n1p==0)n1p=i;
						else
						{
							data[i].id=0;
							dataleft--;
						}
					}
				data[n1p].num=0;
			}
		}
	}
}
/*----------------------------------------------------------------------------------------------------------*/
void tranfunvar(const int x,const int y)
{
	int i;
	for(i=x+1;i<y;i++)
		if(data[i].id>0&&data[i].cla>=9&&data[i].cla<=12)
		{
			data[i].num=var[data[i].cla-8].temp;
			data[i].cla=1;
		}
}
/*----------------------------------------------------------------------------------------------------------*/
void funvarnext()
{
	int funuselist[5]={0},funusesum=0,i,j,t=0;
	for(i=4;i>=1;i--)
		if(var[i].use==1)
		{
			funuselist[funusesum]=i;
			funusesum++;
		}
	for(i=0;i<funusesum;i++)
	{
		if(var[funuselist[i]].temp+var[funuselist[i]].step<=var[funuselist[i]].end)
		{
			var[funuselist[i]].temp+=var[funuselist[i]].step;
			t=1;
			for(j=0;j<i;j++)
				var[funuselist[j]].temp=var[funuselist[j]].start;
			break;
		}
	}
	if(t==0)var[funuselist[0]].temp=var[funuselist[0]].end+(1e+300);
}
/*----------------------------------------------------------------------------------------------------------*/
void fundatacopy(const double f,const int z)
{
	int i;
	for(i=1;i<=4;i++)
		fopti[z].x[i]=var[i].temp;
	fopti[z].f=f;
}
/*----------------------------------------------------------------------------------------------------------*/
bool funend()
{
	int i;
	for(i=1;i<=4;i++)
		if(var[i].use==1&&var[i].temp>var[i].end)return 0;
	return 1;
}
/*----------------------------------------------------------------------------------------------------------*/
bool include(const int x,const int z)
{
	if(z==1)
		if(x==1||x==2)return 1;
	if(z==2)
		if(x==3||x==4||x==6||x==7||x==103)return 1;
	if(z==3)
		if(x==5)return 1;
	if(z==4)
		if(x==45||x==46||x==47||x==95)return 1;
	if(z==5)
		if(x==75||x==76||x==77)return 1;
	if(z==6)
		if(x==79||x==80)return 1;
	if(z==7)
		if(x==29||x==30)return 1;
	if(z==8)
		if(x==121)return 1;
	if(z==9)
		if(x==107||x==108||x==109||x==111||x==112)return 1;
	return 0;
}
/*----------------------------------------------------------------------------------------------------------*/
void countmatrix(const int x,const int y,const int z)
{
	int i,j,k,temp1,temp2;
	double a[2001]={0},temp;
	bool *t;
	if(z==1||z==2||z==75||z==76||z==77)
	{
		if(data[x].r!=data[y].r||data[x].c!=data[y].c)matherror=1;
		if(data[x].fac==0||data[y].fac==0)data[x].fac=data[y].fac=0;
		for(i=1;i<=data[x].r*data[x].c&&i<=data[y].r*data[y].c;i++)
			data[x].mat[i]=count(data[x].mat[i],data[y].mat[i],z,&matherror,degreemode);
		delete []data[y].mat;
		return;
	}
	if(z==3)
	{
		if(data[x].r==1&&data[x].c==2&&data[y].r==1&&data[y].c==2)
		{
			a[1]=data[x].mat[1]*data[y].mat[1]-data[x].mat[2]*data[y].mat[2];
			a[2]=data[x].mat[1]*data[y].mat[2]+data[x].mat[2]*data[y].mat[1];
			data[x].mat[1]=a[1];
			data[x].mat[2]=a[2];
			data[x].fac=data[y].fac=1;
			delete []data[y].mat;
			return;
		}
		if(data[x].r==1&&data[x].c==3&&data[y].r==1&&data[y].c==3)
		{
			a[1]=data[x].mat[2]*data[y].mat[3]-data[x].mat[3]*data[y].mat[2];
			a[2]=data[x].mat[3]*data[y].mat[1]-data[x].mat[1]*data[y].mat[3];
			a[3]=data[x].mat[1]*data[y].mat[2]-data[x].mat[2]*data[y].mat[1];
			for(i=1;i<=3;i++)
				data[x].mat[i]=a[i];
			delete []data[y].mat;
			return;
		}
		if(data[y].r!=data[x].c||data[x].r*data[y].c>2000)
		{
			delete []data[x].mat;
			delete []data[y].mat;
			matherr(x);
			return;
		}
		for(i=1;i<=data[x].r;i++)
			for(j=1;j<=data[y].c;j++)
			{
				temp=0;
				for(k=1;k<=data[x].c;k++)
				{
					temp1=(i-1)*data[x].c+k;
					temp2=(k-1)*data[y].c+j;
					temp+=data[x].mat[temp1]*data[y].mat[temp2];
				}
				a[(i-1)*data[y].c+j]=temp;
			}
		delete []data[x].mat;
		data[x].mat=new double[data[x].r*data[y].c+2];
		for(i=1;i<=data[x].r*data[y].c;i++)
			data[x].mat[i]=a[i];
		data[x].c=data[y].c;
		delete []data[y].mat;
		return;
	}
	if(z==103)
	{
		temp=0;
		if(data[x].r*data[x].c!=data[y].r*data[y].c)matherror=1;
		for(i=1;i<=data[x].r*data[x].c&&i<=data[y].r*data[y].c;i++)
			temp+=data[x].mat[i]*data[y].mat[i];
		delete []data[x].mat;
		delete []data[y].mat;
		data[x].cla=1;
		data[x].num=temp;
		return;
	}
	if(z==4)
	{
		if(data[x].r==1&&data[x].c==2&&data[y].r==1&&data[y].c==2)
		{
			if(data[y].mat[1]==0&&data[y].mat[2]==0)
			{
				delete []data[y].mat;
				delete []data[x].mat;
				matherr(x);
				return;
			}
			a[1]=data[x].mat[1]*data[y].mat[1]+data[x].mat[2]*data[y].mat[2];
			a[2]=data[x].mat[2]*data[y].mat[1]-data[x].mat[1]*data[y].mat[2];
			data[x].fac=data[y].fac=1;
			for(i=1;i<=2;i++)
				data[x].mat[i]=a[i]/(data[y].mat[1]*data[y].mat[1]+data[y].mat[2]*data[y].mat[2]);
			delete []data[y].mat;
			return;
		}
		invmat(y);
		countmatrix(x,y,3);
		return;
	}
	if(z==107)
	{
		if(data[x].cla!=16)
		{
			data[x].cla=16;
			data[x].fac=0;
			if(data[y].cla==16)data[x].r=data[x].c=data[y].r;
			else data[x].r=data[x].c=1;
			data[x].mat=new double[data[x].r*data[x].r+2];
			for(i=0;i<=data[x].r*data[x].r;i++)
				data[x].mat[i]=0;
			for(i=1;i<=data[x].r*data[x].r;i+=data[x].r+1)
				data[x].mat[i]=data[x].num;
		}
		if(data[y].cla!=16)
		{
			data[y].cla=16;
			data[y].fac=0;
			if(data[x].cla==16)data[y].r=data[y].c=data[x].r;
			else data[y].r=data[y].c=1;
			data[y].mat=new double[data[y].r*data[y].r+2];
			for(i=0;i<=data[y].r*data[y].r;i++)
				data[y].mat[i]=0;
			for(i=1;i<=data[y].r*data[y].r;i+=data[y].r+1)
				data[y].mat[i]=data[y].num;
		}
		if(data[x].r!=data[y].r||(data[x].r*(data[x].c+data[y].c)>2000))
		{
			delete []data[x].mat;
			delete []data[y].mat;
			matherr(x);
			return;
		}
		j=k=0;
		for(i=1;i<=data[x].r*(data[x].c+data[y].c);i++)
		{
			if(i%(data[x].c+data[y].c)>0&&i%(data[x].c+data[y].c)<=data[x].c)
			{
				j++;
				a[i]=data[x].mat[j];
			}
			else
			{
				k++;
				a[i]=data[y].mat[k];
			}
		}
		delete []data[x].mat;
		data[x].c+=data[y].c;
		data[x].mat=new double[data[x].r*data[x].c+2];
		for(i=0;i<=data[x].r*data[x].c;i++)
			data[x].mat[i]=a[i];
		delete []data[y].mat;
		return;
	}
	if(z==108||z==109)
	{
		if(data[x].cla!=16)
		{
			matherr(x);
			if(data[y].cla==16)delete []data[y].mat;
			return;
		}
		k=z==108?data[x].r:data[x].c;
		t=new bool[k+2];
		for(i=0;i<=k;i++)
			t[i]=1;
		if(data[y].cla!=16)
		{
			if(int(data[y].num)>0&&int(data[y].num)<=k)t[int(data[y].num)]=0;
			else
			{
				delete []t;
				delete []data[x].mat;
				matherr(x);
				return;
			}
		}
		else
		{
			for(i=1;i<=data[y].r*data[y].c;i++)
			{
				if(int(data[y].mat[i])>0&&int(data[y].mat[i])<=k)t[int(data[y].mat[i])]=0;
				else
				{
					delete []t;
					delete []data[y].mat;
					delete []data[x].mat;
					matherr(x);
					return;
				}
			}
			delete []data[y].mat;
		}
		j=0;
		t[0]=t[k];
		for(i=1;i<=data[x].r*data[x].c;i++)
			if((z==108&&t[(i-1)/data[x].c+1])||(z==109&&t[i%data[x].c]))
			{
				j++;
				a[j]=data[x].mat[i];
			}
		if(j==0)
		{
			delete []data[x].mat;
			delete []t;
			matherr(x);
			return;
		}
		if(z==108)data[x].r=j/data[x].c;
		else data[x].c=j/data[x].r;
		for(i=1;i<=data[x].r*data[x].c;i++)
			data[x].mat[i]=a[i];
		delete []t;
		return;
	}
	if(z==111||z==112)
	{
		k=z==111?data[x].r:data[x].c;
		if(data[y].cla!=16)
		{
			if(int(data[y].num)>0&&int(data[y].num)<=k)
			{
				if(z==111)
				{
					j=0;
					for(i=(int(data[y].num)-1)*data[x].c+1;i<=int(data[y].num)*data[x].c;i++)
					{
						j++;
						a[j]=data[x].mat[i];
					}
					data[x].r=1;
					for(i=1;i<=data[x].c;i++)
						data[x].mat[i]=a[i];
					return;
				}
				else if(z==112)
				{
					for(i=1;i<=data[x].r;i++)
						a[i]=data[x].mat[int(data[y].num)+(i-1)*data[x].c];
					data[x].c=1;
					for(i=1;i<=data[x].r;i++)
						data[x].mat[i]=a[i];
					return;
				}
			}
			else
			{
				delete []data[x].mat;
				matherr(x);
				return;
			}
		}
		else
		{
			if((z==111&&data[y].r*data[y].c*data[x].c>2000)||(z==112&&data[y].r*data[y].c*data[x].r>2000))
			{
				delete []data[x].mat;
				delete []data[y].mat;
				matherr(x);
				return;
			}
			temp1=0;
			for(i=1;i<=data[y].r*data[y].c;i++)
			{
				if(int(data[y].mat[i])<=0||int(data[y].mat[i])>k)
				{
					delete []data[x].mat;
					delete []data[y].mat;
					matherr(x);
					return;
				}
				else
				{
					if(z==111)
					{
						for(j=(int(data[y].mat[i])-1)*data[x].c+1;j<=int(data[y].mat[i])*data[x].c;j++)
						{
							temp1++;
							a[temp1]=data[x].mat[j];
						}
					}
					else
					{
						for(j=1;j<=data[x].r;j++)
						{
							temp1++;
							a[temp1]=data[x].mat[int(data[y].mat[i])+(j-1)*data[x].c];
						}
					}
				}
			}
			if(z==111)data[x].r=data[y].r*data[y].c;
			else if(z==112)
			{
				data[x].c=data[x].r;
				data[x].r=data[y].r*data[y].c;
			}
			delete []data[x].mat;
			data[x].mat=new double[data[x].r*data[x].c+2];
			for(i=1;i<=data[x].r*data[x].c;i++)
				data[x].mat[i]=a[i];
			delete []data[y].mat;
			if(z==112)transpose(x);
			return;
		}
	}
	if(z==5)
	{
		if(data[x].r==1&&data[x].c==2&&data[y].r==1&&data[y].c==2)
		{
			data[x].fac=data[y].fac=1;
			complexpower(x,y);
		}
		else
		{
			delete []data[y].mat;
			delete []data[x].mat;
			matherr(x);
		}
		return;
	}
	if(z==121)
	{
		if(data[x].cla==16&&data[x].r==1&&data[x].c==1)
		{
			data[x].num=data[x].mat[1];
			delete []data[x].mat;
		}
		else if(data[x].cla==16&&data[x].r==1&&data[x].c==2)
		{
			data[x].num=sqrt(data[x].mat[1]*data[x].mat[1]+data[x].mat[2]*data[x].mat[2]);
			delete []data[x].mat;
		}
		else if(data[x].cla==16)
		{
			delete []data[x].mat;
			if(data[y].cla==16)delete []data[y].mat;
			matherr(x);
			return;
		}
		if(data[y].cla==1)
		{
			if(degreemode)data[y].num/=(180/PI);
		}
		else if(data[y].cla==16&&data[y].r==1&&data[y].c==1)
		{
			data[y].num=data[y].mat[1];
			if(degreemode)data[y].num/=(180/PI);
			delete []data[y].mat;
		}
		else if(data[y].cla==16&&data[y].r==1&&data[y].c==2)
		{
			data[y].num=planeangle(data[y].mat[1],data[y].mat[2]);
			delete []data[y].mat;
		}
		else if(data[y].cla==16)
		{
			delete []data[y].mat;
			matherr(x);
			return;
		}
		createcomplex(x,data[x].num*cos(data[y].num),data[x].num*sin(data[y].num));
		return;
	}
	delete []data[x].mat;
	matherr(x);
}
/*----------------------------------------------------------------------------------------------------------*/
void matrixpower(const int x,double y)
{
	int i,j,k,l,m,d,temp1,temp2;
	double a[2001]={0},**b,temp,arg,rc;
	bool t[50]={0};
	if(data[x].r==1&&data[x].c==2)
	{
		rc=sqrt(data[x].mat[1]*data[x].mat[1]+data[x].mat[2]*data[x].mat[2]);
		arg=planeangle(data[x].mat[1],data[x].mat[2]);
		if(arg>PI)arg-=2*PI;
		rc=__pow(rc,y,&matherror);
		arg*=y;
		data[x].mat[1]=rc*cos(arg);
		data[x].mat[2]=rc*sin(arg);
		return;
	}
	if(data[x].c!=data[x].r||!isint(y))
	{
		matherr(x);
		delete []data[x].mat;
		return;
	}
	y=gauss(y);
	temp=det(data[x].r,data[x].mat);
	if(y<0)
	{
		invmat(x);
		y=-y;
	}
	m=abs(int(y));
	d=0;
	while(m>0)
	{
		d++;
		t[d]=(m^1+1)%2==1?1:0;
		m>>=1;
	}
	if(d==1)return;
	if(d==0)
	{
		for(i=1;i<=data[x].r*data[x].r;i++)
			data[x].mat[i]=0;
		for(i=1;i<=data[x].r*data[x].r;i+=data[x].r+1)
			data[x].mat[i]=1;
		return;
	}
	b=new double*[d+2];
	for(i=0;i<=d+1;i++)
		b[i]=new double[data[x].r*data[x].r+2];
	for(i=0;i<=data[x].r*data[x].r+1;i++)
		b[1][i]=data[x].mat[i];
	for(i=2;i<=d;i++)
	{
		for(j=1;j<=data[x].r;j++)
			for(k=1;k<=data[x].r;k++)
			{
				temp=0;
				for(l=1;l<=data[x].r;l++)
				{
					temp1=(j-1)*data[x].r+l;
					temp2=(l-1)*data[x].r+k;
					temp+=b[i-1][temp1]*b[i-1][temp2];
				}
				a[(j-1)*data[x].r+k]=temp;
			}
		for(j=1;j<=data[x].r*data[x].r;j++)
			b[i][j]=a[j];
	}
	for(i=1;i<=data[x].r*data[x].r;i++)
		b[0][i]=b[d][i];
	for(i=1;i<d;i++)
	{
		if(!t[i])continue;
		for(j=1;j<=data[x].r;j++)
			for(k=1;k<=data[x].r;k++)
			{
				temp=0;
				for(l=1;l<=data[x].r;l++)
				{
					temp1=(j-1)*data[x].r+l;
					temp2=(l-1)*data[x].r+k;
					temp+=b[0][temp1]*b[i][temp2];
				}
				a[(j-1)*data[x].r+k]=temp;
			}
		for(j=1;j<=data[x].r*data[x].r;j++)
			b[0][j]=a[j];
	}
	for(i=1;i<=data[x].r*data[x].r;i++)
		data[x].mat[i]=b[0][i];
	for(i=0;i<=d+1;i++)
		delete []b[i];
	delete []b;
}
/*----------------------------------------------------------------------------------------------------------*/
void complexpower(const int x,const int y)
{
	int temppos[2]={x,y},i;
	bool tempme=matherror;
	double a,b,c,d,re,im,r,arg;
	if(tempme)return;
	if(data[x].cla==1&&data[y].cla==1)
	{
		r=__pow(data[x].num,data[y].num,&matherror);
		if(!matherror)
		{
			data[x].num=r;
			return;
		}
		else matherror=0;
	}
	for(i=0;i<2;i++)
		if(data[temppos[i]].cla!=16)
		{
			data[temppos[i]].cla=16;
			data[temppos[i]].r=1;
			data[temppos[i]].c=2;
			data[temppos[i]].mat=new double[4];
			data[temppos[i]].mat[1]=data[temppos[i]].num;
			data[temppos[i]].mat[2]=0;
			data[temppos[i]].fac=1;
		}
	if(data[x].r!=1||data[x].c!=2||data[y].r!=1||data[y].c!=2)
	{
		delete []data[x].mat;
		delete []data[y].mat;
		matherr(x);
		return;
	}
	a=data[x].mat[1];
	b=data[x].mat[2];
	c=data[y].mat[1];
	d=data[y].mat[2];
	if(b==0&&a==0)
	{
		delete []data[y].mat;
		if(d!=0)
		{
			delete []data[x].mat;
		    matherr(x);
		    return;
		}
		else if(c==0)
		{
			data[x].mat[1]=1;
			data[x].mat[2]=0;
			return;
		}
		else
		{
			data[x].mat[1]=data[x].mat[2]=0;
			return;
		}
	}
	re=log(sqrt(a*a+b*b));
	im=planeangle(a,b);
	if(im>PI)im-=2*PI;
	r=exp(re*c-im*d);
	arg=re*d+im*c;
	data[x].mat[1]=r*cos(arg);
	data[x].mat[2]=r*sin(arg);
	delete []data[y].mat;
}	
/*----------------------------------------------------------------------------------------------------------*/
inline int deepen(const double x,int depnum)
{
	depth[depnum].temp=var[1].temp;
    depnum++;
	depth[depnum].start=x;
	depth[depnum].end=var[1].temp;
	var[1].temp=x;
	depth[depnum].step=var[1].step=(depth[depnum].end-depth[depnum].start)/1000;
	solve1stnum=1;
	return depnum;
}
/*----------------------------------------------------------------------------------------------------------*/
inline int jump(int depnum)
{
	if(depnum>=1)
	{
		depnum<10?depnum--:depnum=0;
		var[1].temp=depth[depnum].temp;
		var[1].step=depth[depnum].step;
		solve1stnum=1;
		eqlt=2;
		return depnum;
	}
	else if(depnum==0&&(var[1].temp>var[1].end)&&var[0].use)
	{
		var[1].temp=depth[0].temp=var[1].start;
		solve1stnum=1;
		eqlt=2;
		var[1].step=depth[0].step=(depth[0].end-depth[0].start)/69301;
		var[0].use=0;
		return depnum;
	}
	else return depnum;
}
/*----------------------------------------------------------------------------------------------------------*/
bool mulroot(const int n,const int m,const double x,const double y)
{
	if(fabs(y)>1e-3)return false;
	if(n>0)
	{
		int i;
		bool t=0;
		for(i=1;i<=n;i++)
			if(fabs(solvelist[i]-x)<1e-7||(i<=m&&fabs(solvelist[i]-x)<1e-2))
			{
				t=1;
				break;
			}
		if(t)return false;
	}
	return true;
}
/*----------------------------------------------------------------------------------------------------------*/
void newroot(const double x,const double y)
{
	if(fabs(y)>1e-4||solvenum2>=300)return;
	int i;
	bool t=0;
	for(i=1;i<=solvenum2;i++)
		if(fabs(x-bestsolve[i].x[1])<1e-2)
		{
			t=1;
			if(fabs(y)<fabs(bestsolve[i].x[1]))
			{
				bestsolve[i].x[1]=x;
				bestsolve[i].f=y;
			}
		}
	if(!t)
	{
		solvenum2++;
		bestsolve[solvenum2].x[1]=x;
		bestsolve[solvenum2].f=y;
	}
}
/*----------------------------------------------------------------------------------------------------------*/
int solveout(int n,const double x,const double y)
{
	n++;
	solvelist[n]=x;
	if(n==1)cout<<"方程的解:"<<endl;
	cout<<"X"<<n;
	if(n>=10)cout<<" = ";
	else cout<<"  = ";
	if(!fracoutmode)cout<<setw(17)<<setprecision(10)<<x;
	else fracout(x,1,17,10);
	cout<<"         abs(L-R) = "<<setprecision(8)<<fabs(y)<<endl;
	return n;
}
/*----------------------------------------------------------------------------------------------------------*/
void saveanswer()
{
	int i,j,k,l;
	if(memope==0)
	{
		mem[mempos].fac=memans.fac;
		mem[mempos].r=memans.r;
		mem[mempos].c=memans.c;
		for(j=1;j<=mem[mempos].r*mem[mempos].c;j++)
			mem[mempos].mat[j]=memans.mat[j];
	}
	else if(memope==1||memope==2)
	{
		if(mem[mempos].r==memans.r&&mem[mempos].c==memans.c)
		{
			mem[mempos].fac=maxnum(mem[mempos].fac,memans.fac);
			for(j=1;j<=memans.r*memans.c;j++)
				memope==1?mem[mempos].mat[j]+=memans.mat[j]:mem[mempos].mat[j]-=memans.mat[j];
		}
		else if(mem[mempos].r==1&&memans.r==1&&(mem[mempos].c==1||mem[mempos].c==2)&&(memans.c==1||memans.c==2))
		{
			mem[mempos].fac=mem[mempos].r=1;
			mem[mempos].c=2;
			for(j=1;j<=memans.r*memans.c;j++)
				memope==1?mem[mempos].mat[j]+=memans.mat[j]:mem[mempos].mat[j]-=memans.mat[j];
		}
		else if(!funmode&&!solvemode)cout<<"数学错误！无法赋值！\n";
	}
	else if(memope==3)
	{
		if((memans.r==1&&memans.c==1)||(memans.r==1&&memans.c==2&&memans.fac==1&&fabs(memans.mat[2])<1e-15))
		{
			for(i=1;i<=mem[mempos].r*mem[mempos].c;i++)
				mem[mempos].mat[i]*=memans.mat[1];
		}
		else if((mem[mempos].r==1&&mem[mempos].c==1)||(mem[mempos].r==1&&mem[mempos].c==2&&mem[mempos].fac==1&&fabs(mem[mempos].mat[2])<1e-15))
		{
			mem[mempos].c=2;
			mem[0].mat[1]=mem[mempos].mat[1];
			for(i=1;i<=2;i++)
				mem[mempos].mat[i]=memans.mat[i]*mem[0].mat[1];
		}
		else if(memans.r==1&&memans.c==2&&mem[mempos].r==1&&mem[mempos].c==2)
		{
			mem[mempos].fac=1;
			mem[0].mat[1]=mem[mempos].mat[1]*memans.mat[1]-mem[mempos].mat[2]*memans.mat[2];
			mem[mempos].mat[2]=mem[mempos].mat[2]*memans.mat[1]+mem[mempos].mat[1]*memans.mat[2];
			mem[mempos].mat[1]=mem[0].mat[1];
		}
		else if(mem[mempos].c==memans.r&&(mem[mempos].r*memans.c<=2000))
		{
			double temp;
			for(i=1;i<=mem[mempos].r;i++)
				for(j=1;j<=memans.c;j++)
				{
					temp=0;
					for(k=1;k<=memans.c;k++)
						temp+=mem[mempos].mat[(i-1)*memans.r+k]*memans.mat[(k-1)*memans.c+j];
					mem[0].mat[(i-1)*memans.c+j]=temp;
				}
			mem[mempos].c=memans.c;
			for(i=1;i<=mem[mempos].r*mem[mempos].c;i++)
				mem[mempos].mat[i]=mem[0].mat[i];
		}
		else if(!funmode&&!solvemode)cout<<"数学错误！无法赋值！\n";
	}
	else if(memope==4)
	{
		if(memans.r==mem[mempos].r&&(memans.r*(mem[mempos].c+memans.c)<=2000))
		{
			j=k=0;
			mem[mempos].fac=2;
    		for(l=1;l<=memans.r*(mem[mempos].c+memans.c);l++)
			{
				if(l%(mem[mempos].c+memans.c)>0&&l%(mem[mempos].c+memans.c)<=mem[mempos].c)
				{
					j++;
					mem[0].mat[l]=mem[mempos].mat[j];
				}
				else
				{
					k++;
					mem[0].mat[l]=memans.mat[k];
				}
			}
			mem[mempos].c+=memans.c;
			for(j=1;j<=mem[mempos].r*mem[mempos].c;j++)
				mem[mempos].mat[j]=mem[0].mat[j];
		}
		else if(!funmode&&!solvemode)cout<<"数学错误！无法赋值！\n";
	}
}
/*----------------------------------------------------------------------------------------------------------*/
void writememory()
{
	ofstream outfile;
	outfile.open("CountFever.sav",ios::out);
	int i,j;
	outfile<<memans.fac<<endl<<memans.r<<endl<<memans.c<<endl;
	for(i=1;i<=memans.r*memans.c&&i<=2000;i++)
		outfile<<memans.mat[i]<<endl;
	for(i=1;i<=26;i++)
	{
		outfile<<mem[i].fac<<endl<<mem[i].r<<endl<<mem[i].c<<endl;
		for(j=1;j<=mem[i].r*mem[i].c&&j<=2000;j++)
			outfile<<mem[i].mat[j]<<endl;
	}
	for(i=1;i<=6;i++)
	{
		for(j=0;j<strlen(ordinfun[i]);j++)
			outfile<<ordinfun[i][j];
		outfile<<endl;
	}
	for(i=0;i<=9;i++)
	{
		for(j=0;j<strlen(memstr[i]);j++)
			outfile<<memstr[i][j];
		outfile<<endl;
	}
	outfile.close();
}
/*----------------------------------------------------------------------------------------------------------*/
void debugoutput()
{
	int i;
	if(funmode==1)cout<<"funmode on"<<endl;
	for(i=1;i<=datasum;i++)
	{
		if(data[i].id!=0)
		{
			if(data[i].cla==1)cout<<' '<<data[i].num<<' ';
			else if(data[i].cla==8)cout<<" k ";
			else if(data[i].cla==9)cout<<" X ";
			else if(data[i].cla==10)cout<<" Y ";
			else if(data[i].cla==11)cout<<" Z ";
			else if(data[i].cla==12)cout<<" W ";
			else if(data[i].cla==16)cout<<" M("<<data[i].r<<"*"<<data[i].c<<") ";
			else if(data[i].ope<=100)
			{
				if(data[i].ope==1)cout<<" + ";
				else if(data[i].ope==2)cout<<" - ";
				else if(data[i].ope==3)cout<<" * ";
				else if(data[i].ope==4)cout<<" / ";
				else if(data[i].ope==5)cout<<" ^ ";
				else if(data[i].ope==6)cout<<" div ";
				else if(data[i].ope==7)cout<<" mod ";
				else if(data[i].ope==8)cout<<" ! ";
				else if(data[i].ope==9)cout<<"+";
				else if(data[i].ope==10)cout<<"-";
				else if(data[i].ope==11)cout<<" ( ";
				else if(data[i].ope==12)cout<<" ln( ";
				else if(data[i].ope==13)cout<<" abs( ";
				else if(data[i].ope==14)cout<<" ) ";
				else if(data[i].ope==15)cout<<" sqrt( ";
				else if(data[i].ope==16)cout<<" sin( ";
				else if(data[i].ope==17)cout<<" cos( ";
				else if(data[i].ope==18)cout<<" tan( ";
				else if(data[i].ope==19)cout<<" lg( ";
				else if(data[i].ope==20)cout<<" exp( ";
				else if(data[i].ope==21)cout<<" ctg( ";
				else if(data[i].ope==22)cout<<" asin( ";
				else if(data[i].ope==23)cout<<" acos( ";
				else if(data[i].ope==24)cout<<" atan( ";
				else if(data[i].ope==25)cout<<" sh( ";
				else if(data[i].ope==26)cout<<" ch( ";
				else if(data[i].ope==27)cout<<" th( ";
				else if(data[i].ope==28)cout<<" log( ";
				else if(data[i].ope==29)cout<<" C ";
				else if(data[i].ope==30)cout<<" P ";
				else if(data[i].ope==31)cout<<" gcd( ";
				else if(data[i].ope==32)cout<<" lcm( ";
				else if(data[i].ope==33)cout<<" , ";
				else if(data[i].ope==34)cout<<" [ ";
				else if(data[i].ope==35)cout<<" { ";
				else if(data[i].ope==36)cout<<" ran( ";
				else if(data[i].ope==37)cout<<" sum(k= ";
				else if(data[i].ope==38)cout<<" pro(k= ";
				else if(data[i].ope==39)cout<<" I( ";
				else if(data[i].ope==40)cout<<" % ";
				else if(data[i].ope==41)cout<<" ^2 ";
				else if(data[i].ope==42)cout<<" ^3 ";
				else if(data[i].ope==43)cout<<" max( ";
				else if(data[i].ope==44)cout<<" min( ";
				else if(data[i].ope==45)cout<<" & ";
				else if(data[i].ope==46)cout<<" | ";
				else if(data[i].ope==47)cout<<" xor ";
				else if(data[i].ope==48)cout<<" _ ";
				else if(data[i].ope==49)cout<<" < ";
				else if(data[i].ope==50)cout<<" > ";
				else if(data[i].ope==51)cout<<" == ";
				else if(data[i].ope==52)cout<<" <= ";
				else if(data[i].ope==53)cout<<" >= ";
				else if(data[i].ope==54)cout<<" ≠ ";
				else if(data[i].ope==55)cout<<" sec( ";
				else if(data[i].ope==56)cout<<" csc( ";
				else if(data[i].ope==57)cout<<" ash( ";
				else if(data[i].ope==58)cout<<" ach( ";
				else if(data[i].ope==59)cout<<" ath( ";
				else if(data[i].ope==60)cout<<" ? ";
				else if(data[i].ope==61)cout<<" dX) ";
				else if(data[i].ope==62)cout<<" dY) ";
				else if(data[i].ope==63)cout<<" dZ) ";
				else if(data[i].ope==64)cout<<" dW) ";
				else if(data[i].ope==65)cout<<" f( ";
				else if(data[i].ope==66)cout<<" g( ";
				else if(data[i].ope==67)cout<<" h( ";
				else if(data[i].ope==68)cout<<" F( ";
				else if(data[i].ope==69)cout<<" G( ";
				else if(data[i].ope==70)cout<<" H( ";
				else if(data[i].ope==71)cout<<" fib( ";
				else if(data[i].ope==72)cout<<" sgn( ";
				else if(data[i].ope==73)cout<<" root( ";
				else if(data[i].ope==74)cout<<" dis( ";
				else if(data[i].ope==75)cout<<" b& ";
				else if(data[i].ope==76)cout<<" b| ";
				else if(data[i].ope==77)cout<<" b^ ";
				else if(data[i].ope==78)cout<<" ~ ";
				else if(data[i].ope==79)cout<<" << ";
				else if(data[i].ope==80)cout<<" >> ";
				else if(data[i].ope==81)cout<<" rnd( ";
				else if(data[i].ope==82)cout<<" gpa( ";
				else if(data[i].ope==83)cout<<" actg( ";
				else if(data[i].ope==84)cout<<" sum( ";
				else if(data[i].ope==85)cout<<" sum2( ";
				else if(data[i].ope==86)cout<<" on( ";
				else if(data[i].ope==87)cout<<" on-1( ";
				else if(data[i].ope==88)cout<<" avg( ";
				else if(data[i].ope==89)cout<<" mid( ";
				else if(data[i].ope==90)cout<<" n( ";
				else if(data[i].ope==91)cout<<" avg2( ";
				else if(data[i].ope==92)cout<<" ga( ";
				else if(data[i].ope==93)cout<<" ha( ";
				else if(data[i].ope==94)cout<<" pro( ";
				else if(data[i].ope==95)cout<<" xnor ";
				else if(data[i].ope==96)cout<<" diff( ";
				else if(data[i].ope==97)cout<<" det( ";
				else if(data[i].ope==98)cout<<" M( ";
				else if(data[i].ope==99)cout<<" ; ";
				else if(data[i].ope==100)cout<<" : ";
			}
			else if(data[i].ope<=200)
			{
				if(data[i].ope==101)cout<<" ' ";
				else if(data[i].ope==102)cout<<" "<<char(34)<<" ";
				else if(data[i].ope==103)cout<<" ` ";
				else if(data[i].ope==104)cout<<" rank( ";
				else if(data[i].ope==105)cout<<" diag( ";
				else if(data[i].ope==106)cout<<" tr( ";
				else if(data[i].ope==107)cout<<" $ ";
				else if(data[i].ope==108)cout<<" dr ";
				else if(data[i].ope==109)cout<<" dc ";
				else if(data[i].ope==110)cout<<" size( ";
				else if(data[i].ope==111)cout<<" gr ";
				else if(data[i].ope==112)cout<<" gc ";
				else if(data[i].ope==113)cout<<" a( ";
				else if(data[i].ope==114)cout<<" re( ";
				else if(data[i].ope==115)cout<<" im( ";
				else if(data[i].ope==116)cout<<" sumabs( ";
				else if(data[i].ope==117)cout<<" solve( ";
				else if(data[i].ope==118)cout<<" floor( ";
				else if(data[i].ope==119)cout<<" ceil( ";
				else if(data[i].ope==120)cout<<" c( ";
				else if(data[i].ope==121)cout<<" \\ ";
				else if(data[i].ope==122)cout<<" r( ";
				else if(data[i].ope==123)cout<<" up( ";
				else if(data[i].ope==124)cout<<" dowm( ";
				else if(data[i].ope==126)cout<<" cod( ";
				else if(data[i].ope==127)cout<<" sumxy"<<data[i].fac<<"( ";
				else if(data[i].ope==128)cout<<" fc"<<data[i].fac<<"( ";
				else if(data[i].ope==129)cout<<" A[ ";
			}
		}
	}
	cout<<endl;
}
/*----------------------------------------------------------------------------------------------------------*/
int main()
{
	int i,j=0,k,wtime,ftime,f1stnum,temppos=0,length,depnum,stime,solvenum,solvenum3,outfill;
	bool solved,alwaysmatherr,cmpable,cmpable2,complexfracout;
	double eqlx,eps;
	const string dayofweek[7]={"星期日","星期一","星期二","星期三","星期四","星期五","星期六"};
	time_t calctime;
	initializeandreadmemory();
	randominitialize();
	cout<<setprecision(12);
	while(1)
	{
		var[1].temp=var[2].temp=var[3].temp=var[4].temp=fopti[0].f=fopti[1].f=ftime=stime=depnum=solvenum=solvenum2=solvenum3=syntaxerror=solved=cmpable2=plotting=0;
		solve1stnum=alwaysmatherr=cmpable=1;
		eqlt=2;
		for(i=1;i<=4;i++)
			fopti[0].x[i]=fopti[1].x[i]=0;
		do
		{
		    inputexpression();
		}while(specialmode);
		if(quit)break;
		calctime=clock();
		translate();
		if(funmode==1&&solvemode==0)
		{
			cout<<"   N |";
			j=f1stnum=0;
			for(i=1;i<=4;i++)
				if(var[i].use==1)
				{
					cout<<"               "<<char(i<=3?i+87:87)<<" |";
					j++;
				}
			cout<<" ";
			temppos=0;
			length=strlen(input);
			if(ordinfunmode==1)
			    for(i=0;i<length;i++)
				    if(input[i]=='?')
					{
					    temppos=i+1;
					    break;
					}
			for(i=0;i<2*j+2;i++)
				cout<<input[i+temppos];
			cout<<endl<<"-----|";
			for(i=1;i<=j;i++)
				cout<<"-----------------|";
			cout<<"--------------------------------------------\n";
		}
		do
		{
			randominitialize();
			eps=solvemode?1e-13:1e-11;
			if(funmode==1)
				for(i=1;i<=4;i++)
					if(var[i].use==1)
						if(fabs(var[i].temp-gauss(var[i].temp))<eps)var[i].temp=gauss(var[i].temp);
		    dataleft=datasum;
		    wtime=matherror=0;
			ondisplay=1;
		    if(syntaxerror)cout<<"表达式语法错误！\n";
		    else
			{
		        wtime=maincalc();
				if(funmode==1&&ondisplay==1&&solvemode==0)
				{
					if(!plotting)cout<<setw(4)<<ftime+1<<" |";
				    for(i=1;i<=4;i++)
					    if(var[i].use==1)
						{
							if(plotting&&funplot)plotdata[ftime].x=var[i].temp;
							else if(!fracoutmode)cout<<setw(16)<<setprecision(8)<<var[i].temp<<" |";
							else
							{
								fracout(var[i].temp,1,16,8);
								cout<<" |";
							}
						}
				}
				if((wtime>=400||syntaxerror)&&ondisplay==1)funmode?cout<<" 表达式语法错误！\n":cout<<"表达式语法错误！\n";
			    else if(matherror==1&&ondisplay==1)
				{
					if(solvemode)
					{
						if(solve1stnum)
						{
							eqlx=var[1].temp;
							if(eqlt==0)stime++;
							else stime=0;
							eqlt=0;
						}
						else if(eqlt!=0)
						{
							depnum=deepen(eqlx,depnum);
							stime=0;
						}
					}
					else if(plotting)plotdata[ftime].re=plotdata[ftime].im=plotdata[ftime].t=0;
					else funmode?cout<<" 数学错误！\n":cout<<"数学错误！\n";
				}
		        else if(ondisplay==1)
				{
				    if(datasum==0&&mempos>0&&mempos<27)
					{
					    datasum=1;
					    data[1].id=1;
						if(memans.fac==0)
						{
						    data[1].cla=1;
					        data[1].num=memans.mat[1];
						}
						else
						{
							if(memans.fac==0)data[datasum].num=memans.mat[1];
					        else
							{
						        if(memans.fac==1)data[datasum].fac=1;
						        else data[datasum].fac=2;
						        data[datasum].cla=16;
						        data[datasum].r=memans.r;
						        data[datasum].c=memans.c;
						        data[datasum].mat=new double[memans.r*memans.c+2];
						        for(j=1;j<=memans.r*memans.c;j++)
							    data[datasum].mat[j]=memans.mat[j];
							}
						}
					}
					syntaxerror=1;
		            for(i=1;i<=datasum;i++)
			            if(data[i].id>0&&(data[i].cla==16||data[i].cla==1||data[i].cla==17))
						{
							if(data[i].cla==16&&data[i].c==1&&data[i].r==1)oneomatrix(i);
							if(plotting)
							{
								syntaxerror=0;
								if(data[i].cla==17||(data[i].cla==16&&data[i].fac!=1))
								{
									plotdata[ftime].re=plotdata[ftime].im=plotdata[ftime].t=0;
									if(data[i].cla==16)delete []data[i].mat;
								}
								else if(data[i].cla==16)
								{
									if(_isnan(data[i].mat[1])!=0||(1.0e+306)/data[i].mat[1]==0||_isnan(data[i].mat[2])!=0||(1.0e+306)/data[i].mat[2]==0)plotdata[ftime].re=plotdata[ftime].im=plotdata[ftime].t=0;
									else
									{
									    plotdata[ftime].t=1;
									    plotdata[ftime].re=data[i].mat[1];
									    plotdata[ftime].im=data[i].mat[2];
									}
									delete []data[i].mat;
								}
								else if(_isnan(data[i].num)!=0||(1.0e+306)/data[i].num==0)plotdata[ftime].re=plotdata[ftime].im=plotdata[ftime].t=0;
								else
								{
									plotdata[ftime].t=1;
									plotdata[ftime].re=data[i].num;
									plotdata[ftime].im=0;
								}
								break;
							}
							if(solvemode&&data[i].cla==17)
							{
								data[i].num=1;
								matherror=1;
							}
							if(solvemode&&data[i].cla==16)
							{
								data[i].cla=1;
								data[i].num=0;
								for(j=1;j<=data[i].r*data[i].c;j++)
									data[i].num+=data[i].mat[j]*data[i].mat[j];
								data[i].num=sqrt(data[i].num);
								delete []data[i].mat;
							}
							if(!funmode&&!solvemode&&data[i].cla==16)
							{
								syntaxerror=0;
								for(j=1;j<=data[i].r*data[i].c;j++)
									if(_isnan(data[i].mat[j])!=0||(1.0e+306)/data[i].mat[j]==0)
									{
										matherror=1;
										break;
									}
								if(matherror)
								{
									cout<<"数学错误！\n";
									delete []data[i].mat;
									break;
								}
								if(data[i].c>70)
								{
									cout<<"矩阵太大无法输出！\n";
									delete []data[i].mat;
									break;
								}
								cout<<"答案:             (计算耗时"<<clock()-calctime+1<<"毫秒)\n";
								if(data[i].r==1&&data[i].c==2&&data[i].fac==1)
								{
									complexfracout=complexout(data[i].mat[1],data[i].mat[2]);
									cout<<endl;
									if(complexfracout)
									{
										complexdecout(data[i].mat[1],data[i].mat[2]);
										cout<<endl;
									}
								}
								for(j=1;j<=data[i].r;j++)
								{
									if(data[i].r==1&&data[i].c==2&&data[i].fac==1)break;
									if(data[i].r==1)cout<<"（";
									else if(j==1)cout<<"┌";
									else if(j==data[i].r)cout<<"└";
									else cout<<"│";
									for(k=1;k<=data[i].c;k++)
									{
										if(fabs(data[i].mat[(j-1)*data[i].c+k])<1e-14)data[i].mat[(j-1)*data[i].c+k]=0;
										if(data[i].c>65)outfill=14;
										else if(data[i].c>61)outfill=15;
										else if(data[i].c>54)outfill=16;
										else outfill=17;
										if(fracoutmode)
										{
											fracout(data[i].mat[(j-1)*data[i].c+k],1,outfill,outfill==17?10:outfill-8);
											if(data[i].c<=54)cout<<" ";
										}
										else
										{
										    if(data[i].c>65)cout<<setprecision(6)<<setw(14)<<data[i].mat[(j-1)*data[i].c+k];
										    else if(data[i].c>61)cout<<setprecision(7)<<setw(15)<<data[i].mat[(j-1)*data[i].c+k];
										    else if(data[i].c>54)cout<<setprecision(8)<<setw(16)<<data[i].mat[(j-1)*data[i].c+k];
										    else cout<<setprecision(10)<<setw(17)<<data[i].mat[(j-1)*data[i].c+k]<<" ";
										}
										outfill=0;
									}
									if(data[i].r==1)cout<<"    ）";
									else if(j==1)cout<<"    ┐";
									else if(j==data[i].r)cout<<"    ┘";
									else cout<<"    │";
									if(j==data[i].r)cout<<data[i].r<<"×"<<data[i].c<<endl;
									else cout<<endl;
								}
								memans.r=data[i].r;
								memans.c=data[i].c;
								memans.fac=data[i].fac==1?1:2;
								for(j=1;j<=memans.r*memans.c&&j<=2000;j++)
									memans.mat[j]=data[i].mat[j];
								if(mempos>0&&mempos<27&&!solvemode)saveanswer();
								delete []data[i].mat;
								break;
							}
							if(!solvemode&&!funmode&&data[i].cla==17)
							{
								data[i].num=countdate(&data[i].fac,&data[i].r,&data[i].c,0,&matherror);
								cout<<"答案:             (计算耗时"<<clock()-calctime+1<<"毫秒)\n";
								cout<<data[i].fac<<"年"<<data[i].r<<"月"<<data[i].c<<"日 "<<dayofweek[int(data[i].num)]<<endl;
								syntaxerror=0;
								memans.fac=0;
								memans.r=memans.c=1;
								memans.mat[1]=data[i].num;
								data[i].cla=1;
								break;
							}
					        if(_isnan(data[i].num)!=0||(1.0e+306)/data[i].num==0)
							{
								if(solvemode)
								{
						            if(solve1stnum)
									{
							            eqlx=var[1].temp;
										if(eqlt==0)stime++;
										else stime=0;
					            		eqlt=0;
									}
						            else if(eqlt!=0)
									{
							            depnum=deepen(eqlx,depnum);
										stime=0;
									}
								}
					            else funmode?cout<<" 数学错误！\n":cout<<"数学错误！\n";
							}
						    else 
							{
								if(funmode==1&&fabs(data[i].num)<1e-11&&!solvemode)data[i].num=0;
								if(funmode==1)
								{
									alwaysmatherr=0;
									if(solvemode)newroot(var[1].temp,data[i].num);
									if(solvemode==0)
									{
										if(data[i].cla==16)
										{
											memans.fac=data[i].fac==1?1:2;
											memans.r=data[i].r;
											memans.c=data[i].c;
											for(j=1;j<=data[i].r*data[i].c&&j<=2000;j++)
												memans.mat[j]=data[i].mat[j];
										}
										if(data[i].cla==16&&data[i].r*data[i].c>50)
										{
											delete []data[i].mat;
											data[i].num=0;
											cmpable=0;
											cout<<" 矩阵太大无法输出！\n";
										}
										else if(data[i].cla==16&&data[i].r==1&&data[i].c==2&&data[i].fac==1)
										{
											if(fabs(data[i].mat[2])<1e-11)
											{
												cmpable2=1;
												if(fabs(data[i].mat[1])<1e-11)data[i].mat[1]=0;
												data[i].num=data[i].mat[1];
												if(fracoutmode)
												{
													fracout(data[i].mat[1],1,18,10);
													cout<<endl;
												}
												else cout<<setw(18)<<setprecision(10)<<data[i].mat[1]<<endl;
											}
											else if(fabs(data[i].mat[1])<1e-11)
											{
												data[i].num=cmpable=0;
												if(fabs(fabs(data[i].mat[2])-1)>1e-9)
												{
													if(data[i].mat[2]>0)cout<<" ";
													if(fracoutmode)
													{
														cout<<"                  ";
														fracout(data[i].mat[2],2,NULL,10);
														cout<<"i\n";
													}
													else cout<<"                  "<<setprecision(10)<<data[i].mat[2]<<"i"<<endl;
												}
												else if(data[i].mat[2]>0)cout<<"                   i\n";
												else cout<<"                  -i\n";
											}
											else
											{
												data[i].num=cmpable=0;
												if(fracoutmode)fracout(data[i].mat[1],1,18,10);
												else cout<<setw(18)<<setprecision(10)<<data[i].mat[1];
												if(data[i].mat[2]>0)cout<<"+";
												else cout<<"-";
												data[i].mat[2]=fabs(data[i].mat[2]);
												if(fabs(data[i].mat[2]-1)>1e-9)
												{
													if(fracoutmode)
													{
														fracout(data[i].mat[2],2,NULL,10);
														cout<<"i\n";
													}
													else cout<<setprecision(10)<<data[i].mat[2]<<"i\n";
												}
												else cout<<"i\n";
											}
											delete []data[i].mat;
										}
										else if(data[i].cla==16)
										{
											cmpable=0;
											cout<<" [ ";
											for(j=1;j<=data[i].r*data[i].c;j++)
											{
												if(fabs(data[i].mat[j])<1e-11)data[i].mat[j]=0;
												if(fracoutmode)fracout(data[i].mat[j],1,15,7);
												else cout<<setprecision(7)<<setw(15)<<data[i].mat[j];
												if(j==data[i].r*data[i].c)cout<<" ] "<<data[i].r<<"×"<<data[i].c<<endl;
												else if(j%data[i].c==0)cout<<" ; ";
												else cout<<" , ";
											}
											delete []data[i].mat;
											data[i].num=0;
										}
										else if(data[i].cla==17)
										{
											cmpable=0;
											data[i].num=countdate(&data[i].fac,&data[i].r,&data[i].c,0,&matherror);
											memans.fac=0;
								            memans.r=memans.c=1;
							            	memans.mat[1]=data[i].num;
											data[i].cla=1;
											cout<<setw(8)<<data[i].fac<<" 年 "<<setw(2)<<data[i].r<<" 月 "<<setw(2)<<data[i].c<<" 日  "<<dayofweek[int(data[i].num)]<<endl;
										}
										else 
										{
											cmpable2=1;
											if(fracoutmode)
											{
												fracout(data[i].num,1,18,10);
												cout<<endl;
											}
											else cout<<setprecision(10)<<setw(18)<<data[i].num<<endl;
										}
									}
									else if(fabs(data[i].num)<1e-15)
									{
										if(mulroot(solvenum,solvenum,var[1].temp,data[i].num))
										{
											solvenum=solveout(solvenum,var[1].temp,data[i].num);
										    solved=1;
										}
										depnum=jump(depnum);
									}
									else if(solve1stnum==1&&eqlt==0)
									{
										depnum=deepen(eqlx,depnum);
										stime=0;
									}
									else if(solve1stnum==1)
									{
										if(data[i].num>0)eqlt=1;
										else eqlt=-1;
										eqlx=var[1].temp;
										solve1stnum=0;
										stime=0;
									}
									else if(eqlt*data[i].num<0)
									{
										depnum=deepen(eqlx,depnum);
										stime=0;
									}
									else
									{
										eqlx=var[1].temp;
										solve1stnum=0;
										stime++;
									}
									if(f1stnum==0)
									{
										f1stnum=1;
										fundatacopy(data[i].num,0);
										fundatacopy(data[i].num,1);
									}
									else
									{
										if(data[i].num>fopti[0].f)fundatacopy(data[i].num,0);
										if(data[i].num<fopti[1].f)fundatacopy(data[i].num,1);
									}
								}
								else
								{
									cout<<"答案:             (计算耗时"<<clock()-calctime+1<<"毫秒)\n";
									fracout(data[i].num,4,NULL,12);
									cout<<setprecision(12)<<data[i].num<<endl;
								}
							    if(data[i].cla==1)
								{
									memans.r=memans.c=1;
									memans.fac=0;
									memans.mat[1]=data[i].num;
								}
							    if(mempos>0&&mempos<27&&!solvemode)saveanswer();
							}
							syntaxerror=0;
				            break;
						}
					if(syntaxerror)funmode?cout<<" 表达式语法错误！\n":cout<<"表达式语法错误！\n";
				}
			}
			if((syntaxerror||wtime>=400)&&solvemode)break;
			if(solvemode)
			{
				var[1].temp+=var[1].step;
				if((var[1].temp>depth[depnum].end)||(depnum>1&&stime>1500)||depnum>=10)
				{
					if(depnum>1&&stime>1500&&memans.mat[1]<1e-8&&mulroot(solvenum,solvenum,var[1].temp,memans.mat[1]))
					{
						solvenum=solveout(solvenum,var[1].temp,memans.mat[1]);
						solved=1;
					}
					depnum=jump(depnum);
					stime=0;
				}
			}
			else if(funmode)funvarnext();
			if(solvemode&&solvenum==solvesum)
			{
				memans.mat[1]=var[1].temp;
				break;
			}
			if(funmode==1)translate();
			if(ondisplay==1)ftime++;
			if(funmode&&!solvemode&&funplot&&!plotting&&(!funend()||ftime>=9000))
			{
				for(j=1;j<=4;j++)
					if(var[j].use==1)
					{
						var[j].end-=1e-9;
						var[j].temp=var[j].start;
						var[j].step=(var[j].end-var[j].start)/640;
						var[j].end+=1e-9;
						break;
					}
				ftime=0;
				for(j=0;j<=640;j++)
					plotdata[j].x=plotdata[j].re=plotdata[j].im=plotdata[j].t=0;
				plotting=1;
			}
		}while(funmode==1&&funend()&&((ftime<9000&&solvemode==0)||(ftime<350000&&solvemode)));
		if(solvemode&&!syntaxerror)
		{
			solvenum3=solvenum;
			for(i=1;i<=solvenum2&&solvenum<solvesum;i++)
			{
				if(mulroot(solvenum,solvenum3,bestsolve[i].x[1],bestsolve[i].f))
				{
					solvenum=solveout(solvenum,bestsolve[i].x[1],bestsolve[i].f);
					solved=1;
				}
			}
			if(alwaysmatherr)cout<<"数学错误！\n";
			else if(!solved)cout<<"方程很可能无解！\n";
		}
		if(solvemode||funmode)
		{
			memans.r=memans.c=1;
			memans.mat[1]=memans.fac=0;
			solvesum=solvenum;
		}
		writememory();
		if(ftime==9000)cout<<"数据太多！ \n";
		if(funmode==1&&solvemode==0)
		{
			if(cmpable&&cmpable2)
			{
			    cout<<"以上函数值中的最大值: ";
			    for(i=1;i<=4;i++)
				    if(var[i].use==1)
					{
					    cout<<"  "<<char(i<=3?i+87:87)<<"= ";
						if(fracoutmode)fracout(fopti[0].x[i],2,NULL,10);
					    else cout<<setprecision(10)<<fopti[0].x[i];
					}
			    cout<<"   "<<input[temppos]<<"= ";
				if(fracoutmode)fracout(fopti[0].f,2,NULL,10);
				else cout<<setprecision(10)<<fopti[0].f;
			    cout<<"\n以上函数值中的最小值: ";
			    for(i=1;i<=4;i++)
				    if(var[i].use==1)
					{
					    cout<<"  "<<char(i<=3?i+87:87)<<"= ";
					    if(fracoutmode)fracout(fopti[1].x[i],2,NULL,10);
					    else cout<<setprecision(10)<<fopti[1].x[i];
					}
			    cout<<"   "<<input[temppos]<<"= ";
				if(fracoutmode)fracout(fopti[1].f,2,NULL,10);
				else cout<<setprecision(10)<<fopti[1].f;
				cout<<endl;
			}
			else cout<<"无法比较以上函数值的大小！\n";
		}
		if(plotting)plot(input,ordinoy,codoy);
	}
	cout<<"正在储存数据...\n\n";
	Sleep(250);
	return 0;
}