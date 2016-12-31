#include <graphics.h>
#include <conio.h>
using namespace std;
int plot(char *title,const bool ordinoy,const double codoy); //绘制一元函数图象
struct countfeverfunctionvaluedataforplot
{
	double x,re,im,rediff,imdiff,rei,imi;
	int re0,im0,rediff0,imdiff0,rei0,imi0;
	bool t,diff;
}plotdata[641];
/*----------------------------------------------------------------------------------------------------------*/
int plot(char *title,const bool ordinoy,const double codoy)
{
	const char title1[100]=" - Count Fever 1.08 中文版";
	double ox,oy,step,value[641]={0},maxy,miny,resum,imsum;
	int i,j,total=0,codx,cody;
	step=plotdata[1].x-plotdata[0].x;
	for(i=0;i<=640;i++)
	{
		if(plotdata[i].t)
		{
			value[total]=plotdata[i].re;
			total++;
		}
		if(i>0&&i<800&&plotdata[i].t&&plotdata[i-1].t&&plotdata[i+1].t)
		{
			plotdata[i].diff=1;
			plotdata[i].rediff=(plotdata[i+1].re-plotdata[i-1].re)/(2*step);
			plotdata[i].imdiff=(plotdata[i+1].im-plotdata[i-1].im)/(2*step);
		}
		else plotdata[i].diff=0;
	}
	for(i=0;i<total;i++)
		for(j=i+1;j<total;j++)
			if(value[i]>value[j])
			{
				oy=value[i];
				value[i]=value[j];
				value[j]=oy;
			}
	if(total==0)return -1;
	const COLORREF COLORWHITE=RGB(255,255,255);
	const COLORREF COLORRED=RGB(255,0,0);
	const COLORREF COLORGREEN=RGB(0,255,0);
	const COLORREF COLORBLUE=RGB(20,20,255);
	const COLORREF COLORGREY=RGB(63,63,63);
	const COLORREF COLORYELLOW=RGB(255,127,0);
	const COLORREF COLORLIGHTBLUE=RGB(31,191,159);
	const COLORREF COLORPURPLE=RGB(255,127,127);
	const COLORREF COLORDARKGREEN=RGB(31,95,0);
	const COLORREF COLORLIGHTGREY=RGB(191,191,191);
	ox=(plotdata[0].x+plotdata[640].x)/2;
	if(!ordinoy)
	{
	    if(value[total/5]>=-320*step&&value[total-total/5-1]<=320*step)oy=0;
	    else oy=__int64((value[total/2]+value[(total-1)/2])/2);
	}else oy=codoy;
	maxy=oy+320*step;
	miny=oy-320*step;
	plotdata[320].rei=plotdata[320].imi=0;
	resum=imsum=0;
	for(i=321;i<=640;i++)
		if(plotdata[i].t)
		{
			resum+=plotdata[i].re;
			imsum+=plotdata[i].im;
			plotdata[i].rei=resum*step;
			plotdata[i].imi=imsum*step;
		}
	resum=imsum=0;
	for(i=319;i>=0;i--)
		if(plotdata[i].t)
		{
			resum+=plotdata[i].re;
			imsum+=plotdata[i].im;
			plotdata[i].rei=-resum*step;
			plotdata[i].imi=-imsum*step;
		}
	initgraph(800,640,SHOWCONSOLE);
	HWND hWnd=GetHWnd();
	SetWindowText(hWnd,strcat(title,title1));
	cleardevice();
	setcolor(COLORGREY);
	setlinestyle(PS_DASH,1);
	for(i=0;i<640;i+=64)
	{
		if(i==320)continue;
		line(i,0,i,640);
	}
	for(i=0;i<640;i+=64)
	{
		if(i==320)continue;
		line(0,i,640,i);
	}
	setcolor(COLORLIGHTGREY);
	if(plotdata[0].x<0&&plotdata[640].x>0&&miny<0&&maxy>0)
	{
		codx=int(640*plotdata[0].x/(plotdata[0].x-plotdata[640].x));
		cody=int(640*maxy/(maxy-miny));
		if(codx!=320||cody!=320)
		{
		    line(0,cody,640,cody);
		    line(codx,0,codx,640);
			RECT r0={codx-15,cody+15,codx-2,cody+2};
			setcolor(COLORWHITE);
			drawtext("O",&r0,DT_RIGHT|DT_VCENTER|DT_SINGLELINE);
		}
	}
	else if(plotdata[0].x<0&&plotdata[640].x>0)
	{
		codx=int(640*plotdata[0].x/(plotdata[0].x-plotdata[640].x));
		line(codx,0,codx,640);
	}
	else if(miny<0&&maxy>0)
	{
		cody=int(640*maxy/(maxy-miny));
		line(0,cody,640,cody);
	}
	setlinestyle(PS_SOLID,1);
	setcolor(COLORWHITE);
	line(0,320,640,320);
	line(635,325,640,320);
	line(635,315,640,320);
	line(320,0,320,640);
	line(320,0,315,5);
	line(320,0,325,5);
	for(i=0;i<=640;i++)
	{
		if(!plotdata[i].t)continue;
		if(320-(plotdata[i].imi-oy)/step>2147481647)plotdata[i].imi0=2147481647;
		else if(320-(plotdata[i].imi-oy)/step<-2147481648)plotdata[i].imi0=-2147481648;
		else plotdata[i].imi0=320-(plotdata[i].imi-oy)/step;
		setcolor(COLORDARKGREEN);
		if(i==0||(i>0&&plotdata[i-1].t==0)||fabs(plotdata[i].imi-plotdata[i-1].imi)/(maxy-miny)>3)putpixel(i,plotdata[i].imi0,COLORGREEN);
		else line(i-1,plotdata[i-1].imi0,i,plotdata[i].imi0);
		if(320-(plotdata[i].rei-oy)/step>2147481647)plotdata[i].rei0=2147481647;
		else if(320-(plotdata[i].rei-oy)/step<-2147481648)plotdata[i].rei0=-2147481648;
		else plotdata[i].rei0=320-(plotdata[i].rei-oy)/step;
		setcolor(COLORPURPLE);
		if(i==0||(i>0&&plotdata[i-1].t==0)||fabs(plotdata[i].rei-plotdata[i-1].rei)/(maxy-miny)>3)putpixel(i,plotdata[i].rei0,COLORRED);
		else line(i-1,plotdata[i-1].rei0,i,plotdata[i].rei0);
		if(plotdata[i].diff)
		{
			if(320-(plotdata[i].imdiff-oy)/step>2147481647)plotdata[i].imdiff0=2147481647;
			else if(320-(plotdata[i].imdiff-oy)/step<-2147481648)plotdata[i].imdiff0=-2147481648;
			else plotdata[i].imdiff0=320-(plotdata[i].imdiff-oy)/step;
			setcolor(COLORLIGHTBLUE);
			if(i==0||(i>0&&plotdata[i-1].diff==0)||fabs(plotdata[i].imdiff-plotdata[i-1].imdiff)/(maxy-miny)>0.5)putpixel(i,plotdata[i].imdiff0,COLORLIGHTBLUE);
			else line(i-1,plotdata[i-1].imdiff0,i,plotdata[i].imdiff0);
			if(320-(plotdata[i].rediff-oy)/step>2147481647)plotdata[i].rediff0=2147481647;
			else if(320-(plotdata[i].rediff-oy)/step<-2147481648)plotdata[i].rediff0=-2147481648;
			else plotdata[i].rediff0=320-(plotdata[i].rediff-oy)/step;
			setcolor(COLORYELLOW);
			if(i==0||(i>0&&plotdata[i-1].diff==0)||fabs(plotdata[i].rediff-plotdata[i-1].rediff)/(maxy-miny)>0.5)putpixel(i,plotdata[i].rediff0,COLORYELLOW);
			else line(i-1,plotdata[i-1].rediff0,i,plotdata[i].rediff0);
		}
		if(320-(plotdata[i].im-oy)/step>2147481647)plotdata[i].im0=2147481647;
		else if(320-(plotdata[i].im-oy)/step<-2147481648)plotdata[i].im0=-2147481648;
		else plotdata[i].im0=320-(plotdata[i].im-oy)/step;
		setcolor(COLORGREEN);
		if(i==0||(i>0&&plotdata[i-1].t==0)||fabs(plotdata[i].im-plotdata[i-1].im)/(maxy-miny)>1.5)putpixel(i,plotdata[i].im0,COLORGREEN);
		else line(i-1,plotdata[i-1].im0,i,plotdata[i].im0);
		if(320-(plotdata[i].re-oy)/step>2147481647)plotdata[i].re0=2147481647;
		else if(320-(plotdata[i].re-oy)/step<-2147481648)plotdata[i].re0=-2147481648;
		else plotdata[i].re0=320-(plotdata[i].re-oy)/step;
		setcolor(COLORRED);
		if(i==0||(i>0&&plotdata[i-1].t==0)||fabs(plotdata[i].re-plotdata[i-1].re)/(maxy-miny)>1.5)putpixel(i,plotdata[i].re0,COLORRED);
		else line(i-1,plotdata[i-1].re0,i,plotdata[i].re0);
	}
	ofstream outfile;
	outfile.open("temp.txt",ios::out);
	outfile<<setprecision(6)<<plotdata[0].x<<endl<<plotdata[640].x<<endl<<miny<<endl<<maxy<<endl<<"("<<ox<<","<<oy<<")";
	outfile.close();
	ifstream infile;
	infile.open("temp.txt",ios::in);
	char tempch1[20],tempch2[20],tempch3[20],tempch4[20],tempch5[20];
	infile>>tempch1>>tempch2>>tempch3>>tempch4>>tempch5;
	infile.close();
	remove("temp.txt");
	RECT r1={3,340,90,320};
	RECT r2={550,340,635,320};
	RECT r3={230,640,315,620};
	RECT r4={230,20,315,0};
	RECT r5={170,340,318,320};
	RECT r6={3,640,100,620};
	setbkmode(OPAQUE|TRANSPARENT);
	setcolor(COLORWHITE);
	drawtext(tempch1,&r1,DT_LEFT|DT_VCENTER|DT_SINGLELINE);
	drawtext(tempch2,&r2,DT_RIGHT|DT_VCENTER|DT_SINGLELINE);
	drawtext(tempch3,&r3,DT_RIGHT|DT_VCENTER|DT_SINGLELINE);
	drawtext(tempch4,&r4,DT_RIGHT|DT_VCENTER|DT_SINGLELINE);
	drawtext(tempch5,&r5,DT_RIGHT|DT_VCENTER|DT_SINGLELINE);
	setcolor(COLORWHITE);
	drawtext("按任意键返回",&r6,DT_LEFT|DT_VCENTER|DT_SINGLELINE);
//	getch();
//	closegraph();
	return 0;
}
