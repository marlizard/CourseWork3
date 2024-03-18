#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

typedef double real;

//Номер используемой функции
#define function 4

//Используется ли генератор тестов и сетки
#define generator 1

//Используется ли шум
#define noise 0

//Величина шума в процентах
#define magnitude 1

class Spline {
public:
	//Количество интервалов сплайна
	int n;

	//Координаты узлов сплайна
	vector<real> knot;

	//Количество наборов
	int m;

	//Координаты точек наборов
	vector<real> x;

	//Значения в точках набора
	vector<real> f;

	//Значения сплайна в точках
	vector<real> y;

	//Глобальная матрица
	vector<vector<real>> A;

	//Вектор правой части
	vector<real> b;

	//Вектор коэффициентов сглаживания
	vector<real> q;

	//Вектор весов для регуляции близости сплайна
	vector<real> omega;

	//Параметр регуляризации
	real beta = 10e-18;

	//Искомая функция
	real func(real x);

	//Первая эрмитова базисная функция
	real first_hermit(real x, int i);

	//Вторая эрмитова базисная функция
	real second_hermit(real x, int i);

	//Третья эрмитова базисная функция
	real third_hermit(real x, int i);

	//Четвертая эрмитова базисная функция
	real forth_hermit(real x, int i);

	//Генератор сетки и тестов
	void grid_and_test_generator();

	//Ввод данных
	void input();

	//Посторение портрета матрицы
	void matrix_portrait();

	//Учет локальных матриц
	void local_matrix(int i);

	//Сборка глобальной матрицы
	void global_matrix();

	//Решаетль методом Гаусса
	void Gauss_method();

	//Расчет значений сплайна
	void spline_values();

	//Вывод результата
	void output();
};


real add_noise(real volume) {
	real start = volume / 100 * (100 - magnitude);
	real end = volume / 100 * (100 + magnitude);

	//return start + rand() / volume / 10000 * (end - start);
	return start + volume / rand() / 100 * (end - start);
}


real Spline::func(real x) {
	switch (function) {
	case 0: {
		return 5;
		break;
	}
	case 1: {
		return x;
		break;
	}
	case 2: {
		return pow(x, 2);
		break;
	}
	case 3: {
		return pow(x, 3);
		break;
	}
	case 4: {
		return pow(x, 4);
		break;
	}
	case 5: {
		return cos(x) + 2 * pow(sin(x), 2);
		break;
	}
	}
}


real Spline::first_hermit(real x, int i) {
	real ksi = (x - knot[i]) / (knot[i + 1] - knot[i]);

	return 1 - 3 * pow(ksi, 2) + 2 * pow(ksi, 3);
}


real Spline::second_hermit(real x, int i) {
	real h = knot[i + 1] - knot[i];
	real ksi = (x - knot[i]) / h;

	return h * (ksi - 2 * pow(ksi, 2) + pow(ksi, 3));
}


real Spline::third_hermit(real x, int i) {
	real ksi = (x - knot[i]) / (knot[i + 1] - knot[i]);

	return 3 * pow(ksi, 2) - 2 * pow(ksi, 3);
}


real Spline::forth_hermit(real x, int i) {
	real h = knot[i + 1] - knot[i];
	real ksi = (x - knot[i]) / h;

	return h * (- pow(ksi, 2) + pow(ksi, 3));
}


void Spline::grid_and_test_generator() {
	ifstream generate("generate.txt");
	int finite_el, set, start, finish;
	generate >> finite_el >> set;
	generate >> start >> finish;
	generate.close();

	real h_fe = finish - start;
	h_fe /= finite_el;
	real knots = start;

	ofstream grid("grid.txt");
	grid << finite_el << "\n";
	for (int i = 0; i <= finite_el; i++) {
		grid << knots << "  ";
		knots += h_fe;
	}
	grid.close();

	real h_s = finish - start;
	h_s /= (set - 1);
	real points = start;

	ofstream values("values.txt");
	values << set << "\n";
	if (noise == 1) {
		for (int i = 0; i < set; i++) {
			values << fixed;
			values << setprecision(14) << points << "  " << add_noise(func(points)) << "  " << 1 << "\n";
			points += h_s;
		}
	}
	else {
		for (int i = 0; i < set; i++) {
			values << fixed;
			values << setprecision(16) << points << "  " << func(points) << "  " << 1 << "\n";
			points += h_s;
		}
	}
	values.close();
}


void Spline::input() {
	ifstream values("values.txt");
	values >> m;
	x.resize(m);
	f.resize(m);
	omega.resize(m);
	y.resize(m, 0);

	for (int i = 0; i < m; i++) {
		values >> x[i] >> f[i] >> omega[i];
	}

	values.close();

	ifstream grid("grid.txt");
	grid >> n;
	knot.resize(n + 1);
	for (int i = 0; i <= n; i++) grid >> knot[i];
	grid.close();
}


void Spline::matrix_portrait() {
	A.resize(2 * n + 2, vector<real>(2 * n + 2, 0));

	b.resize(2 * n + 2, 0);
	q.resize(2 * n + 2, 0);
}


void Spline::local_matrix(int i) {
	real h = knot[i + 1] - knot[i];
	
	//Прибавка третьего компонента
	A[2 * i][2 * i] += beta * 12 / pow(h, 3);
	A[2 * i + 1][2 * i + 1] += beta * 4 / h;
	A[2 * i + 2][2 * i + 2] += beta * 12 / pow(h, 3);
	A[2 * i + 3][2 * i + 3] += beta * 4 / h;

	A[2 * i + 1][2 * i] += beta * 6 / pow(h, 2);
	A[2 * i + 2][2 * i] -= beta * 12 / pow(h, 3);
	A[2 * i + 2][2 * i + 1] -= beta * 6 / pow(h, 2);
	A[2 * i + 3][2 * i] += beta * 6 / pow(h, 2);
	A[2 * i + 3][2 * i + 1] += beta * 2 / h;
	A[2 * i + 3][2 * i + 2] -= beta * 6 / pow(h, 2);

	A[2 * i][2 * i + 1] += beta * 6 / pow(h, 2);
	A[2 * i][2 * i + 2] -= beta * 12 / pow(h, 3);
	A[2 * i + 1][2 * i + 2] -= beta * 6 / pow(h, 2);
	A[2 * i][2 * i + 3] += beta * 6 / pow(h, 2);
	A[2 * i + 1][2 * i + 3] += beta * 2 / h;
	A[2 * i + 2][2 * i + 3] -= beta * 6 / pow(h, 2);

	for (int k = 0; k < m; k++) {
		if (x[k] >= knot[i]) {
			if (x[k] < knot[i + 1]) {

				real first = first_hermit(x[k], i);
				real second = second_hermit(x[k], i);
				real third = third_hermit(x[k], i);
				real forth = forth_hermit(x[k], i);
				real w = omega[k];

				//Прибавка первого компонента
				A[2 * i][2 * i] += w * pow(first, 2);
				A[2 * i + 1][2 * i + 1] += w * pow(second, 2);
				A[2 * i + 2][2 * i + 2] += w * pow(third, 2);
				A[2 * i + 3][2 * i + 3] += w * pow(forth, 2);

				A[2 * i + 1][2 * i] += w * first * second;
				A[2 * i + 2][2 * i] += w * first * third;
				A[2 * i + 2][2 * i + 1] += w * second * third;
				A[2 * i + 3][2 * i] += w * first * forth; 
				A[2 * i + 3][2 * i + 1] += w * second * forth; 
				A[2 * i + 3][2 * i + 2] += w * third * forth; 

				A[2 * i][2 * i + 1] += w * first * second;
				A[2 * i][2 * i + 2] += w * first * third;
				A[2 * i + 1][2 * i + 2] += w * second * third;
				A[2 * i][2 * i + 3] += w * first * forth;
				A[2 * i + 1][2 * i + 3] += w * second * forth;
				A[2 * i + 2][2 * i + 3] += w * third * forth;

				//Вектор правой части
				b[2 * i] += w * f[k] * first;
				b[2 * i + 1] += w * f[k] * second;
				b[2 * i + 2] += w * f[k] * third;
				b[2 * i + 3] += w * f[k] * forth;
			}
			else break;
		}
	}

	if (i == knot.size() - 2 && x.back() == knot.back()) {

		int k = x.size() - 1;

		real first = first_hermit(x[k], i);
		real second = second_hermit(x[k], i);
		real third = third_hermit(x[k], i);
		real forth = forth_hermit(x[k], i);
		real w = omega[k];

		A[2 * i][2 * i] += w * pow(first, 2);
		A[2 * i + 1][2 * i + 1] += w * pow(second, 2);
		A[2 * i + 2][2 * i + 2] += w * pow(third, 2);
		A[2 * i + 3][2 * i + 3] += w * pow(forth, 2);

		A[2 * i + 1][2 * i] += w * first * second;
		A[2 * i + 2][2 * i] += w * first * third;
		A[2 * i + 2][2 * i + 1] += w * second * third;
		A[2 * i + 3][2 * i] += w * first * forth;
		A[2 * i + 3][2 * i + 1] += w * second * forth;
		A[2 * i + 3][2 * i + 2] += w * third * forth;

		A[2 * i][2 * i + 1] += w * first * second;
		A[2 * i][2 * i + 2] += w * first * third;
		A[2 * i + 1][2 * i + 2] += w * second * third;
		A[2 * i][2 * i + 3] += w * first * forth;
		A[2 * i + 1][2 * i + 3] += w * second * forth;
		A[2 * i + 2][2 * i + 3] += w * third * forth;

		//Вектор правой части
		b[2 * i] += w * f[k] * first;
		b[2 * i + 1] += w * f[k] * second;
		b[2 * i + 2] += w * f[k] * third;
		b[2 * i + 3] += w * f[k] * forth;
	}
}


void Spline::global_matrix() {
	for (int i = 0; i < n; i++) local_matrix(i);
}


void Spline::Gauss_method() {
	for (int i = 1; i < 2 * n + 2; i++)
		for (int j = i; j < 2 * n + 2; j++) {
			real s = A[j][i - 1] / A[i - 1][i - 1];
			for (int k = 0; k < 2 * n + 2; k++)
				A[j][k] = A[j][k] - s * A[i - 1][k];
			b[j] = b[j] - s * b[i - 1];
		}

	for (int i = 2 * n + 1; i >= 0; i--) {
		real s = 0;
		for (int j = i + 1; j < 2 * n + 2; j++)
			s += A[i][j] * b[j];
		b[i] -= s;
		b[i] /= A[i][i];
	}

	for (int i = 0; i < q.size(); i++) q[i] = b[i];
}


void Spline::spline_values() {
	for (int i = 0; i < m; i++) {
		for (int k = 0; k < n; k++) {
			if (x[i] >= knot[k] && x[i] <= knot[k + 1])
				y[i] = q[2 * k] * first_hermit(x[i], k) +
					   q[2 * k + 1] * second_hermit(x[i], k) +
					   q[2 * k + 2] * third_hermit(x[i], k) +
					   q[2 * k + 3] * forth_hermit(x[i], k);
		}
	}
}


void Spline::output() {
	ofstream out("output.txt");
	out << "\t\tx\t\t\t\t\tu\t\t\t\t\ty\t\t\t\t|u - y|\t\t\t\t\tf\t\t\t\t|u - f|\n";
	for (int i = 0; i < m; i++) {
		out << fixed << setprecision(10);
		out <<  scientific << x[i] << "\t" << func(x[i]) << "\t" 
			<< y[i] << "\t" << abs(y[i] - func(x[i])) << "\t" 
			<< f[i] << "\t" << abs(f[i] - func(x[i])) << "\n";
	}

	out.close();
}


int main() {
	Spline m;
	if(generator) m.grid_and_test_generator();
	m.input();
	m.matrix_portrait();
	m.global_matrix();
	m.Gauss_method();
	m.spline_values();
	m.output();
}