//Программа решает дифуравнение методом Рунге-Кутты 4го порядка, 
//подбирая для расчётов шаг заданной точности методом двойного пересчёта,
//после решается интеграл методом Симпсона, под интегралом - y^2

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <conio.h>
#include <vector>
using namespace std;
double f(double, double);// правая часть диф.уравнения
double DoubleCounting(double, double, double, double, double, vector<double>&, vector<double>&);// функция двойного пересчёта
double RungeKuttaRank4(double, double, double);//вычисляет n+1 значение y с заданным n-м и заданным шагом методом Рунге-Кутты 4го порядка
double yx(double, double, double, vector<double>&, vector<double>&);// y(x) на основе значений полученный методом Рунге-Кутты и интерполяции промежуточных значений
double Simpson(double start,double end,double, double, vector<double>&, vector<double>&);//вычисление интеграла методом Симпсона

double Simpson(//вычисление интеграла методом Симпсона
		double start,		//начало интервала интегрирования
		double end,			//конец интервала интегрирования
		double h, 			//шаг метода Симпсона
		double hRK,			//итоговый шаг метода РГ(Рунге-Кутты)
		vector<double> &vxn,	//значения xi
		vector<double> &vyn){	//значения yi
	double Int; //значение интеграла
	Int = vyn[0]*vyn[0] + vyn[vyn.size()-1]*vyn[vyn.size()-1]; //начальное значение интегральной суммы - под интегралом у нас y^2
	int n = (end-start)/h;
	double xn = start;//текущее значение
	double yt = vyn[0];
	for(int i=1; i<n; ++i){
		xn = xn + h;
		yt = yx(xn, start, hRK,vxn, vyn);
		if(i%2 == 1) 	Int = Int + 4*yt*yt;
		else 			Int = Int + 2*yt*yt;
	}
	cout<<endl;
	Int = Int *h/3;
	return Int;
}

double yx(// y(x) на основе значений полученный методом Рунге-Кутты и интерполяции промежуточных значений
		double x, 				//аргумент искомого значения функции
		double x0, 				//начальный x - начало промежутка, левый край
		double h, 				//итоговое значение шага метода Рунге-Кутты, на котором была достигнута требуемая точность
		vector<double>&vxn, 	//значения xi
		vector<double>&vyn){	//значение yi
	int i = int((x-x0)/h);		//вычисляем примерный номер искомого значения функции
	if(vxn[i]==x) return vyn[i];
	else {
		double y = vyn[i]+(vyn[i+1]-vyn[i])*(x-vxn[i])/h;
		cout<<"Интерполяция в точке x= "<<x<<" : y= "<<y<<endl;
		return y; //если необходимое значение аргумента функции не найдено в таблице значений полученных решением диф.уравнения,
		//применяем линейную интерполяцию
	}
}
double f(double x, double y) // правая часть диф.уравнения
{
	return 1+3*y*sin(x)-y;
}

double RungeKuttaRank4 ( double xn, double yn, double h){ //вычисляет n+1 значение y с заданным n-м и заданным шагом методом Рунге-Кутты 4го порядка
	double k1,k2,k3,k4;
	double h2 = h/2;
	k1=f(xn, yn);
	k2=f(xn+h2, yn+k1*h2);
	k3=f(xn+h2, yn+k2*h2);
	k4=f(xn+h, yn+h/k3);
	xn+=h;
	yn+=h*(k1+4*k3+k4)/6;
	return yn;
}

double DoubleCounting(	// функция двойного пересчёта
						double start, 	//начало отрезка
						double end, 	//конец отрезка
						double y0, 		//начальное значение y
						double h, 		//шаг
						double eps,		//точность
						vector<double> &vxn,	//массив значений x
						vector<double> &vyn		//массив значений y
						){
//	vxn.erase(vxn.begin(),vxn.end());
//	vyn.erase(vxn.begin(),vxn.end());
	double xn = start; // начальный x
	double yn = y0;// начальное значение y
	int n = (end-start)/h;
	if(vxn.empty() && vyn.empty()){
		vxn.push_back(xn);
		vyn.push_back(yn);
		for(int i = 0;i < n; ++i){
			yn = RungeKuttaRank4(xn, yn, h);
			xn = xn + h;
			vxn.push_back(xn);
			vyn.push_back(yn);
		}
	}

	double ynh2 = y0; //пересчитанное значение y(n+1) с половинным шагом
	double xnh2 = start;
	vector<double> vynh2;
	vector<double> vxnh2;
	vxnh2.push_back(xnh2);
	vynh2.push_back(ynh2);
	for(int i = 0;i < 2*n; ++i){
		ynh2 = RungeKuttaRank4(xnh2, ynh2, h/2);
		xnh2 = xnh2 + h/2;
		vxnh2.push_back(xnh2);
		vynh2.push_back(ynh2);
	}
	for(int i = 1; i < n; ++i){
		if( fabs(vyn[i] - vynh2[i*2]) > 15* eps) {//если необходимая точность не достигнута, перезапускаем
			vxn.erase(vxn.begin(), vxn.end());
			vyn.erase(vyn.begin(), vyn.end());
			vxn = vxnh2;
			vyn = vynh2;
			vxnh2.erase(vxnh2.begin(), vxnh2.end());
			vynh2.erase(vynh2.begin(), vynh2.end());
			//if(n< 20)
				return DoubleCounting(start, end, y0, h/2, eps, vxn, vyn);
		}
	}
	return h;
}

int main(void){
	vector<double> vxi; //значения x
	vector<double> vyi; //значения y
	double eps = 0.0001;	//точность для решения диф уравнения
	double y0 = 0.2;		//значение y из начального условия
	double start = 0;		//начало интервала
	double end = 1;			//конец интервала
	double h0 = 0.5;		//стартовый шаг для двойного пересчёта
	double h = DoubleCounting(start, end, y0, h0, eps, vxi, vyi); //итоговое значение шага метода Рунге-Кутты 4го ранга, при котором достигается требуемая точность
	int n = vxi.size();
	cout<<"Решение диф уравнения. Значеия функции y в точках x:"<<endl;
	for(int i=0; i < n; ++i)
		cout<<"x= "<<vxi[i]<<"\ty= "<<vyi[i]<<endl;
	cout<<"h="<<h<<endl<<"Q="<<Simpson(start, end, 0.1, h, vxi, vyi);
	getch();
}
