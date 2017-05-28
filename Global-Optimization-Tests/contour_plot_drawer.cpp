
#include "contour_plot_drawer.h"


void ContourPlotDrawer::draw()
{
	/*
	double xray[50], yray[50], zmat[50][50];
	int n = 50, i, j;
	double  fpi = 3.14159 / 180.0, step, x, y;
	double  zlev;

	step = 360.0 / (n - 1);

	for (i = 0; i < n; i++)
	{
		xray[i] = i * step;
		yray[i] = i * step;
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			x = xray[i] * fpi;
			y = yray[j] * fpi;
			zmat[i][j] = 2 * sin(x) * sin(y);
		}
	}

	Dislin dislin;

	dislin.scrmod("revers");
	dislin.setpag("da4p");
	dislin.metafl("cons");

	dislin.disini();
	dislin.complx();
	dislin.pagera();

	dislin.titlin("Contour Plot", 1);
	dislin.titlin("F(X,Y) = 2 * SIN(X) * SIN(Y)", 3);

	dislin.name("X-axis", "x");
	dislin.name("Y-axis", "y");

	dislin.intax();
	dislin.axspos(450, 2670);
	dislin.graf(0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0);

	dislin.height(30);
	for (i = 0; i < 9; i++)
	{
		zlev = -2.0 + i * 0.5;
		dislin.setclr((i + 1) * 25);
		if (i == 4)
			dislin.labels("none", "contur");
		else
			dislin.labels("float", "contur");

		dislin.contur(xray, n, yray, n, (double *)zmat, zlev);
	}
	dislin.height(50);
	dislin.color("fore");
	dislin.title();

	dislin.disfin();
	*/
}
