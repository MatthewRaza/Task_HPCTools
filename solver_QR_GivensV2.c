#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double *generate_matrix(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);  //Vector column pointer.
    time_t t;
   
    /* Intializes random number generator */
    srand((unsigned) time(&t));
   
    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
};

void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;
    printf("matrix: %s \n", name);

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
};



//transposition of a matrix
void transp(double *matrix,int size)              
{
	double temp;
	int i,j;
	
	for(i=0;i<size;i++)
	{
		for(j=0;j<i;j++)
		{
			temp=matrix[i*size+j]; //Because as it is defined, A is not a matrix n*n but a vector with n² rows
			matrix[i*size+j]=matrix[j*size+i];
			matrix[j*size+i]=temp;
		}
	}
};

//multiplication of 2 matrix
double *mult(double *A,double *B,int size)
{
	int i,j,k;
	double *C = (double *)malloc(sizeof(double) * size * size); //initialisaiton of the matrix multiplication
	
	for(i=0;i<size*size;i++)
	{
		C[i]=0.0;
	}
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			for(k=0;k<size;k++)
			{
				C[i*size+j]=C[i*size+j]+A[i*size+k]*B[k*size+j];
			}
		}
	}
	return C;
};

// Rotation of the matrix
double *mat_rot(int size, int i, int j, double c, double s) //Use on QR_Givens
{
	int k;
	double *M = (double *)malloc(sizeof(double) * size * size);
	
	for(k=0;k<size*size;k++)
	{
		M[k]=1;
	}
	M[i*size+i]=c;
	M[j*size+j]=c;
	M[i*size+j]=-s;
	M[j*size+i]=s;
	return M;
}

// QR method with Givens rotation
//-----------------------------------------------------------

//________________________________To have Q and R from Givens methode_______________________

double *R_Givens(double *R, int size)
{
	int i,j;
	double *G;
	
	G=generate_matrix(size);
	
	for(j=0;j<size-1;j++)
	{
		for(i=size-1;i>j;i--)
		{
			int pos1 = i;
			int pos2 = i-1;
			double r;
			
			r=sqrt(R[pos1*size+j]*R[pos1*size+j]+R[(pos2)*size+j]*R[(pos2)*size+j]);
			
			//Creation of a matrix G 
			double *G = generate_matrix(size);
			
			G = mat_rot(size,pos1,(pos2),R[(pos2)*size+j]/r,R[pos1*size+j]/r);
			R = mult(G,R,size);
		}
	}
	
	if(R[(size-1)*size+(size-1)]<0)
	{
		R[(size-1)*size+(size-1)]= R[(size-1)*size+(size-1)]*(-1);
	}
	return R;
	print_matrix("Rint =",R,size);
}

double *Q_Givens(double *R,double *Q, int size)
{
	int i,j;
	double *G = generate_matrix(size);
	
	for(j=0;j<size-1;j++)
	{
		for(i=size-1;i>j;i--)
		{
			int pos1 = i;
			int pos2 = i-1;
			double r;
			
			r=sqrt(R[pos1*size+j]*R[pos1*size+j]+R[(pos2)*size+j]*R[(pos2)*size+j]);
			
			double *G = generate_matrix(size);
			
			G = mat_rot(size,pos1,(pos2),R[(pos2)*size+j]/r,R[pos1*size+j]/r);
			transp(G,size);
			Q = mult(G,R,size);
		}
	}
	
	//Condition to change the co-factor in the matrix.
	if(R[(size-1)*size+(size-1)]<0)
	{
		for(i=0;i<size;i++)
		{
			Q[i*size+2]= Q[i*size+2]*(-1);
		}
	}
	
	return Q;
	print_matrix("Qint =", Q,size);
}

//________________________________________________________________________________

//Resolution of the linear equation AX=B.
double *resolve(double *R,double *B, int size) //fonctionne
{
	double sum;
	int i,j,k;
	double *X = (double *)malloc(sizeof(double) * size * size);
	

	for(i=0;i<size*size;i++)
	{
		X[i]=0.0;
	}
	for(i=size-1;i>=0;i--)
	{
		for(j=0;j<size;j++)
		{
			if (i==size-1)
			{
				if (R[i*size+i]!=0)
				{
					X[i*size+j]=B[i*size+j]/R[i*size+i];
				}
				else
				{
					X[i*size+j]=0;
				}
			}
			else
			{
				sum=0.0;
				for(k=size-1;k>i;k--)
				{
					sum=sum+R[i*size+k]*X[k*size+i];
				}
				if (R[i*size+i]!=0)
				{
					X[i*size+j]=(B[i*size+j]-sum)/R[i*size+j];
				}
				else
				{
					X[i*size+j]=0.0;
				}
			}
		}
	}
	return X;
}

int check_result(double *bref, double *b, int size) {
    int i;
	
    for(i=0;i<size*size;i++) {
        if (bref[i]!=b[i]) return 0;
    }
    return 1;
}

//---------------------------------Main program --------------------------
void main(int argc, char *argv[])
{
   int size = atoi(argv[1]);
   int i;

   double *R, *Q, *NewR, *NewQ;
   double *A, *Aref;
   double *B, *Bref;

    A = generate_matrix(size);
    Aref = generate_matrix(size);        
    B = generate_matrix(size);
    Bref = generate_matrix(size);
        
    print_matrix("A :", A, size);
    print_matrix("B :", B, size);
	
	// Initialization of Q and R
	R=generate_matrix(size);
	Q=generate_matrix(size);
	//---------------------------------
	for(i=0;i<size*size;i++) 
	{
		R[i]=A[i];
		Q[i]=0.0;
	}
	for(i=0;i<size;i++)
	{
		Q[i*size+i]=1.0;
	}
	print_matrix("Q before = ",Q,size);
	print_matrix("R before = ",R,size);
	
	//Change the matrix Q and R by Givens method for resolve AX=B.
	//--------------------------------
	NewR = R_Givens(R,size);
	NewQ = Q_Givens(R,Q,size);
	print_matrix("Q after = ",NewQ,size);
	print_matrix("R after = ",NewR,size);
	transp(Q,size);
	//--------------------------------
	
	double *C;
	C = mult(NewQ,B,size);
	print_matrix(" C : (matricial product QxB)",C,size);
	
	double *X;
	X=resolve(NewR,C,size);
	print_matrix("X, resolution of the linear equation AX=B :",X,size);
}