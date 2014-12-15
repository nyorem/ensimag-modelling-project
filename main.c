#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct s_Point {
	double x, y;
} Point;

#define MAX_VERTICES 10000

#define IDX(i, N) (((i) + (N)) % (N))

Point P1[MAX_VERTICES];
Point P2[MAX_VERTICES];

void PostscriptSetColor( FILE *f, double r, double g, double b)
{
	fprintf( f, "%lf %lf %lf setrgbcolor\n", r, g, b);
}

void PostscriptSetLineWidth( FILE *f, double w)
{
	fprintf( f, "%lf setlinewidth\n", w);
}

void PostscriptHeader( FILE *f)
{
	fprintf( f, "%%!PS\n");	
}

void PostscriptEoF( FILE *f)
{
	fprintf( f, "showpage\n");	
}

void PostscriptOutputClosedPolygon( FILE *f, int N, Point *P)
{
	int i;

	fprintf( f, "%lf %lf moveto\n", P[0].x, P[0].y);
	for (i=1; i<N; i++)
	{
		fprintf( f, "%lf %lf lineto\n", P[i].x, P[i].y);
	}
	fprintf( f, "%lf %lf lineto\n", P[0].x, P[0].y);
	fprintf( f, "stroke\n");
}


void GnuplotOutputClosedPolygon( char *filename, int N, Point *P)
{
	FILE *fout;

	int i;

	fout = fopen( filename, "w");

	for (i=0; i<N; i++)
	{
		fprintf( fout, "%lf %lf\n", P[i].x, P[i].y);
	}
	fprintf( fout, "%lf %lf\n", P[0].x, P[0].y);

	fclose(fout);
}

double AngleThreePoints(  Point *A, Point *B, Point *C)
{
	double V1X, V1Y, V2X, V2Y;
	double N1, N2, ScalarProduct, Angle;

	V1X = B->x - A->x;
	V1Y = B->y - A->y;
	N1 = sqrt( V1X * V1X + V1Y * V1Y);
	V1X /= N1;
	V1Y /= N1;
	

	V2X = C->x - B->x;
	V2Y = C->y - B->y;
	N2 = sqrt( V2X * V2X + V2Y * V2Y);
	V2X /= N2;
	V2Y /= N2;

	ScalarProduct = V1X * V2X + V1Y * V2Y;

	Angle = acos( ScalarProduct);
	if (ScalarProduct < 0.) {Angle = -Angle;}

	Angle = Angle / (N1 + N2);

	return Angle;
}


void ComputeAngle( int N, Point *P, double *Angle)
{
	int i;

	Angle[0] = AngleThreePoints( &P[N-1], &P[0], &P[1]);

	for (i=1; i<N-1; i++)
	{
		Angle[i] = AngleThreePoints( &P[i-1], &P[i], &P[i+1]);
	}

	Angle[i] = AngleThreePoints( &P[N-2], &P[N-1], &P[0]);
}

void GnuplotOutputAngle( char *filename, int N, Point *P)
{
	double Angle[MAX_VERTICES];
	FILE *fout;
	int i;

	ComputeAngle( N, P, Angle);


	fout = fopen( filename, "w");


	for (i=0; i<N; i++)
	{
		fprintf( fout, "%lf %lf\n", ((double) i)/(N-1), Angle[i]);
	}

	fclose(fout);
	
}

void CopyPolygon( int N, Point *P, Point *Q)
{
	int i;

	for (i=0; i<N; i++)
	{
		Q[i] = P[i];
	}
}

double PointDistance( Point *P, Point *Q)
{
	double vx, vy;

	vx = P->x - Q->x;
	vy = P->y - Q->y;

	return(sqrt(vx * vx + vy * vy));
}

double Max( double a, double b)
{
	if (a>b)
	{
		return a;
	}

	return b;
}

// Distance between the "new" points and the previouse ones
double LMaxDistance( int N, Point *P, Point *Q)
{
	double d, di;
	int i;
	Point MidPoint;

	d = 0.;

	for (i=0; i<N-1; i++)
	{
		di = PointDistance( &P[i], &Q[2*i]);
		d = Max( d, di);

		MidPoint.x = 0.5 * (P[i].x + P[i+1].x);
		MidPoint.y = 0.5 * (P[i].y + P[i+1].y);

		di = PointDistance( &MidPoint, &Q[2*i+1]);
		d = Max( d, di);
	}

// CLOSING EDGE

	di = PointDistance( &P[i], &Q[2*i]);
	d = Max( d, di);

	MidPoint.x = 0.5 * (P[i].x + P[0].x);
	MidPoint.y = 0.5 * (P[i].y + P[0].y);

	di = PointDistance( &MidPoint, &Q[2*i+1]);
	d = Max( d, di);

	return(d);
}

void CornerCuttingSubdivisionClosedPolygon( double A, double B, int N, Point *P, Point *Q)
{
	int i;

	for (i=0; i<N-1; i++)
	{
		Q[2*i].x = A * P[i].x + (1. - A) * P[i+1].x;
		Q[2*i].y = A * P[i].y + (1. - A) * P[i+1].y;
		Q[2*i+1].x = B * P[i].x + (1. - B) * P[i+1].x;
		Q[2*i+1].y = B * P[i].y + (1. - B) * P[i+1].y;
	}

	Q[2*i].x = A * P[i].x + (1. - A) * P[0].x;
	Q[2*i].y = A * P[i].y + (1. - A) * P[0].y;
	Q[2*i+1].x = B * P[i].x + (1. - B) * P[0].x;
	Q[2*i+1].y = B * P[i].y + (1. - B) * P[0].y;
}

void ChaikinSubdivisionClosedPolygon( int N, Point *P, Point *Q)
{
	int i;

	for (i=0; i<N-1; i++)
	{
		Q[2*i].x = 3./4. * P[i].x + 1./4. * P[i+1].x;
		Q[2*i].y = 3./4. * P[i].y + 1./4. * P[i+1].y;
		Q[2*i+1].x = 1./4. * P[i].x + 3./4. * P[i+1].x;
		Q[2*i+1].y = 1./4. * P[i].y + 3./4. * P[i+1].y;
	}

	Q[2*i].x = 3./4. * P[i].x + 1./4. * P[0].x;
	Q[2*i].y = 3./4. * P[i].y + 1./4. * P[0].y;
	Q[2*i+1].x = 1./4. * P[i].x + 3./4. * P[0].x;
	Q[2*i+1].y = 1./4. * P[i].y + 3./4. * P[0].y;
}

void GeneralizedFourPointsSchemeSubdivisionClosedPolygon (double epsilon, int N, Point *P, Point *Q)
{
	int i;

	for (i=0; i<N; i++)
	{
		Q[2*i].x = P[i].x;
		Q[2*i].y = P[i].y;
		Q[2*i+1].x = -epsilon * (P[IDX(i-1, N)].x + P[IDX(i+2, N)].x) / 2. + (1 + epsilon) * (P[IDX(i, N)].x + P[IDX(i+1, N)].x) / 2.;
		Q[2*i+1].y = -epsilon * (P[IDX(i-1, N)].y + P[IDX(i+2, N)].y) / 2. + (1 + epsilon) * (P[IDX(i, N)].y + P[IDX(i+1, N)].y) / 2.;
	}
}

void UniformSplinesSubdivisionClosedPolygon (int K, int N, Point *P, Point *Q) {
    // Duplication
    for (int i = 0; i < N; ++i) {
        Q[2*i].x = Q[2*i+1].x = P[i].x;
        Q[2*i].y = Q[2*i+1].y = P[i].y;
    }

    // k averages
    for (int k = 0; k < K; ++k) {
        Point P0;
        P0.x = Q[0].x;
        P0.y = Q[0].y;

        int i;
        for (i = 0; i < 2 * N - 1; ++i) {
            Q[i].x = 0.5 * (Q[i].x + Q[i+1].x);
            Q[i].y = 0.5 * (Q[i].y + Q[i+1].y);
        }

        Q[i].x = 0.5 * (Q[i].x + P0.x);
        Q[i].y = 0.5 * (Q[i].y + P0.y);
    }
}

int main(int argc, char *argv[])
{
	FILE *fout;
	char FileName[100];

    long nbiter = (argc >= 2) ? strtol(argv[1], NULL, 10) : 10;

	int N, i;


	fout = fopen( "PostscriptOutput.ps", "w");

	N = 4;

	P1[0].x = 100.;
	P1[0].y = 300.;
	P1[1].x = 500.;
	P1[1].y = 300.;
	P1[2].x = 500.;
	P1[2].y = 500.;
	P1[3].x = 300.;
	P1[3].y = 700.;

	PostscriptHeader( fout);
	PostscriptSetLineWidth( fout, 1.0);
	PostscriptSetColor( fout, 1., 0., 0.);

	PostscriptOutputClosedPolygon( fout, 4, P1);
	GnuplotOutputClosedPolygon( "Subdivision0.txt", 4, P1);


	for (i=0; i<nbiter; i++)
	{

	/* ChaikinSubdivisionClosedPolygon( N, P1, P2); */
    /* CornerCuttingSubdivisionClosedPolygon( 0.52, 0.41 , N, P1, P2); */
    /* CornerCuttingSubdivisionClosedPolygon( 0.75, 0.25 , N, P1, P2); */
    /* GeneralizedFourPointsSchemeSubdivisionClosedPolygon(1./8., N, P1, P2); */
    UniformSplinesSubdivisionClosedPolygon(3, N, P1, P2);

	PostscriptOutputClosedPolygon( fout, 2 * N, P2);

	sprintf( FileName, "Subdivision%d.txt", i);
	GnuplotOutputClosedPolygon( FileName, 2 * N, P2);
	sprintf( FileName, "AngleSubdivision%d.txt", i);
	GnuplotOutputAngle( FileName, 2 * N, P2);

	printf( "Distance Subdivision: %lf\n", LMaxDistance( N, P1, P2));

	N = N * 2;

	CopyPolygon( N, P2, P1);

	}

	fclose( fout);
}
