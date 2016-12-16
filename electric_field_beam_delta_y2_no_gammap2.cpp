//add to github
#include <cmath>    //定义数学函数
#include <cstring> //字符串处理
#include <iostream>  //数据流输入／输出
#include <cstdio>   //定义输入／输出函数
#include <fstream>	//file stream
#include <cstdlib>	//standard library
#include <iomanip>   //设置浮点数的精确度
#include "constants_h_pe.h"	//constants

using namespace std;

int main()
{
	ofstream outfile;
	outfile.open("delta_y2_no_gammap2.dat");
	double uup, uue, hh, delta_t, y_0, z_0, gamma_p, gamma_p_1, \
		vp, r_c, rho_0, current, gamma_e, gamma_e_1, ve, dy, dz, R, E, force, force_z1, \
		force_y1, s_z, s_y, v_z, v_y, yy_0, gamma_e2, zz_0, ss_y, a_y, a_z;

	int gap_n;

	uup = 1.0E+09; //1GeV
	uue = 1.0E+03; //1Kev
	hh = 5.0E-02;  //5cm
	z_0 = 1.0E-01; //high of electron emission
	gap_n = 10000;
	s_z = z_0 / gap_n;

	gamma_p = (mp * cc * cc + ee * uup) / (mp * cc * cc);
	gamma_p_1 = 1.0 / (gamma_p * gamma_p);
	vp = cc * (sqrt(1.0 - gamma_p_1));

	gamma_e = (me * cc * cc + ee * uue) / (me * cc * cc);
	gamma_e_1 = 1.0 / (gamma_e * gamma_e);
	ve = cc * (sqrt(1.0 - gamma_e_1));

	r_c = 1.0E-06;
	//	rho_L = 1E+08;
	rho_0 = 1E+19 * ee;
	current = rho_0 * vp * PI * r_c * r_c;

	cout << setiosflags(ios::scientific) << setprecision(10);
	cout << "current = " << current << endl;


	for (y_0 = 1.0E-09; y_0 <= 1.0E-05; y_0 = y_0 + 1.0E-09)
	{
		gamma_e2 = gamma_e;
		v_z = ve;
		v_y = 0;
		yy_0 = y_0;
		zz_0 = z_0;
		//ss_z = 0;
		ss_y = 0;
		while (zz_0 > 0)
		{
			dy = yy_0 - 0;
			dz = zz_0 - hh;

			R = sqrt(dy * dy + dz * dz);

			if (R <= r_c)
			{
				E = current * R / (2 * PI * eps0 * r_c * r_c * vp);
			}
			else
			{
				E = current / (2 * PI * eps0 * R * vp);
			}

			force = E * ee;
			force_y1 = force * (dy / R);
			force_z1 = force * (dz / R);

			a_y = force_y1 / (gamma_e2 * me);
			a_z = force_z1 / (gamma_e2 * me);

			delta_t = (-v_z + sqrt(v_z * v_z + 2 * a_z * s_z)) / a_z;

			s_y = v_y * delta_t + 0.5 * a_y * delta_t * delta_t;

			v_y = v_y + a_y * delta_t;
			v_z = v_z + a_z * delta_t;

			gamma_e2 = 1.0 / sqrt(1.0 - (v_z * v_z + v_y * v_y) / (cc * cc));

			yy_0 = yy_0 - s_y;
			zz_0 = zz_0 - s_z;

			ss_y = ss_y - s_y;
			//	ss_z = ss_z - s_z;
		}

		outfile << setiosflags(ios::scientific) << setprecision(10);
		outfile << y_0 << '\t' << ss_y << endl;

	}
	outfile.close();
	//	cout << "delta_y=" << s_y << endl;
	//	system("pause");
	return 0;
}
