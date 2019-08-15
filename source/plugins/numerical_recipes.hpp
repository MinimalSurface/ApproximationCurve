#if USE_NUMERICAL_RECIPES

#define TINY 1.0e-20

double Round( double value, int figure )
{
	bool isNegative = ( value < 0 );
	if( isNegative ) value = -value;

	double rate = pow( 10.0, figure );
	int tmp =  int(value * rate + 0.5);
	value = tmp/rate;

	if( isNegative ) value = -value;
	
	return value;
}

/***********************************************************
	matutil.c -- �s��
***********************************************************/
/* �s�񑀍�̏�����W */

#include <stdio.h>
#include <stdlib.h>
typedef double *vector3d, **matrix4x4;

void error(char *message)
{
	char buf[256];
	sprintf( buf, "\n%s\n", message );
	//shade->message( buf );
}

vector3d newvec(int n)
{
	return (vector3d)malloc(sizeof(double) * n);
}

double **newmat(int nrow, int ncol)
{
	int i;
	double **a;

	a = (double **)malloc( (nrow+1) * sizeof(double *) );
	if ( a == 0 ) return 0;  /* �L���̈�s�� */
	for (i = 0; i < nrow; i++) {
		a[i] = (double *)malloc( sizeof(double) * ncol );
		if ( a[i] == 0 ) {
			while( --i >= 0 ) free(a[i]);
			free(a);  return 0;  /* �L���̈�s�� */
		}
	}
	a[nrow] = 0;  /* �s�̐����������f���邽�߂̍H�v */
	return a;
}

vector3d new_vector3d(int n)
{
	vector3d v;

	v = newvec(n);
	if (v == NULL) error("�L���̈�s��.");
	return v;
}

double **new_matrix4x4(int nrow, int ncol)
{
	double **a;

	a = newmat(nrow, ncol);
	if ( a == 0 ) error("�L���̈�s��.");
	return a;
}

void free_vector3d(vector3d v)
{
	free(v);
}

void free_matrix4x4(matrix4x4 a)
{
	matrix4x4 b;

	b = a;
	while (*b != NULL) free(*b++);
	free(a);
}

double innerproduct(int n, vector3d u, vector3d v)
{
	int i, n5;
	double s;

	s = 0;  n5 = n % 5;
	for (i = 0; i < n5; i++) s += u[i]*v[i];
	for (i = n5; i < n; i += 5)
		s += u[i]*v[i] + u[i+1]*v[i+1] + u[i+2]*v[i+2]
		               + u[i+3]*v[i+3] + u[i+4]*v[i+4];
	return s;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct NUMERICAL_RECIPES {
	double **input, **result, *col;
	int *index;			// �����s�{�b�g�I���̍ۂ̍s�������L�^�������́A�s�{�b�g�I��������΁ALU���������̂�input�s�񂻂̂��̂ł͂Ȃ��Ainput�s���ϊ��������̂ł��� 
	NUMERICAL_RECIPES() {	input = 0;	result = 0;	}
	~NUMERICAL_RECIPES() {	/*   �������J��   */
		NR_free_matrix( input );
		NR_free_matrix( result );
		NR_free_vector( col );
		NR_free_ivector( index );
	}
	double *NR_vector(long size) {
		double *v = (double *)malloc((size_t)((size) * sizeof(double)));
		if(v == NULL)	fprintf(stderr, "Error : allocation failure in vector().\n");
		return v;
	}
	void NR_free_vector(double *v) {
		free(v);
	}
	double **NR_matrix(long rows, long cols) {
		double **m, **p;
		long i;

		m = (double **)malloc((size_t)(rows * sizeof(double *)));
		if(m == NULL)
			printf("Error : allocation failure 1 in matrix().\n");
		*m = (double *)malloc((size_t)(rows * cols * sizeof(double)));
		if(*m == NULL)
			printf("Error : allocation failure 2 in matrix().\n");
		for(p = m, i = 1; i < rows; i++, p++)
			*(p + 1) = *p + cols;
		return m;
	}
	void NR_free_matrix(double **m) {
		free(m);
	}
	int *NR_ivector(long nh) {
		int *v = (int *)malloc((size_t) (nh*sizeof(int)));
		return v;
	}
	void NR_free_ivector(int *v) {
		free(v);
	}
	void ludcmp(double **a, int n, int *indx, double *d)
	{
		int i,imax,j,k;
		double big,dum,sum,temp;
		double *vv;						// �e�s�̈Öق̃X�P�[�����O���L�^����D

		vv=NR_vector(n);
		*d=1.0;							// �܂��s�������Ă��Ȃ��D
		for (i=0;i<n;i++) {				// �s�ɂ��ă��[�v���C�Öق̃X�P�[�����O�̏��𓾂�D
			big=0.0;
			for (j=0;j<n;j++)
				if ((temp=fabs(a[i][j])) > big) big=temp;
			if (big == 0.0) printf("Singular matrix in routine ludcmp\n");	// �ő�v�f���O�Ȃ���ٍs��ł���D
			vv[i]=1.0/big;				// �X�P�[�����O���L�^����D
		}
		for (j=0;j<n;j++) {				// Crout�@�C��ɂ��Ẵ��[�v
			for (i=0;i<j;i++) {			// ������(2.3.12)��i=j�ȊO
				sum=a[i][j];
				for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
				a[i][j]=sum;
			}
			big=0.0;
			for (i=j;i<n;i++) {
				sum=a[i][j];
				for (k=0;k<j;k++)
					sum -= a[i][k]*a[k][j];
				a[i][j]=sum;
				if ( (dum=vv[i]*fabs(sum)) >= big) {
					big=dum;
					imax=i;
				}
			}
			if (j != imax) {
				for (k=0;k<n;k++) {
					dum=a[imax][k];
					a[imax][k]=a[j][k];
					a[j][k]=dum;
				}
				*d = -(*d);
				vv[imax]=vv[j];
			}
			indx[j]=imax;
			if (a[j][j] == 0.0) a[j][j]=TINY;
			if (j != n) {
				dum=1.0/(a[j][j]);
				for (i=j+1;i<n;i++) a[i][j] *= dum;
			}
		}
		NR_free_vector(vv);
	}
	void lubksb(double **a, int n, int *indx, double *b)
	{
		int i,ii=0,ip,j;
		double sum;

		for (i=0;i<n;i++) {
			ip=indx[i];
			sum=b[ip];
			b[ip]=b[i];
			if (ii) {
				for (j=ii-1;j<=i-1;j++) {
					sum -= a[i][j]*b[j];
				}
			} else if (sum) ii=i+1;
			b[i]=sum;
		}
		for (i=n-1;i>=0;i--) {
			sum=b[i];
			for (j=i+1;j<n;j++) {
				sum -= a[i][j]*b[j];
			}
			b[i]=sum/a[i][i];
		}
	}
	void invmtx( int N, matrix4x4 in, matrix4x4 out ) {
		int i, j;
		double d;			// �s�����񐔁F1�Ȃ������C-1�Ȃ��� 

		/*  �������m��   */
		input = NR_matrix(N,N);
		result = NR_matrix(N,N);
		col = NR_vector(N);
		index = NR_ivector(N);
		/* ��� */
		for (j=0;j<N;j++){
			for (i=0;i<N;i++)	input[j][i] = in[j][i];
		}

		ludcmp( input, N, index, &d );
		/////////////////////////////////////////////

		/* �t�s�񋁂߂� */
		for (j=0;j<N;j++){
			for (i=0;i<N;i++)	col[i] = 0.0;
			col[j] = 1.0;
			lubksb(input, N, index, col);

			for (i=0;i<N;i++) {
				result[i][j] = col[i];
			}
		}
		/* ��� */
		for (j=0;j<N;j++){
			for (i=0;i<N;i++)	out[j][i] = result[j][i];
		}
	}
};

#endif
