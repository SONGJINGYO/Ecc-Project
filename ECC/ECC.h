#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<intrin.h>
#include<assert.h>
#include<stdlib.h>
#include<malloc.h>
#include<memory.h>
typedef unsigned int word;
typedef unsigned long long int64_t;

#define Controler 4 // 0: 역원 1:스칼라곱셈 2: 자코비안 3:wNAF(자코비안) 4:Comb
#define ECC_256 8
#define NEGATIVE 1
#define NON_NEGETIVE 0

typedef struct
{
	      // 0 when A >= 0, 1 when A < 0
	word a[ECC_256];
	int sign;
	int carry;

}bigint_st;

typedef struct
{
	      // 0 when A >= 0, 1 when A < 0
	word a[ECC_256+ECC_256];
	int carry;
	int sign;

}bigint_st_mul;

typedef struct {
	bigint_st x;
	bigint_st y;
}ECC_st;

typedef struct {
	bigint_st x;
	bigint_st y;
	bigint_st z;
	int is_infinity;
}ECC_Jacobian_st;

void P_Read(FILE* fp, bigint_st* P);
void P_Read_Mul(FILE* fp, bigint_st_mul* P);
void Show_integer(bigint_st* P);
void Set_integer(bigint_st* P, int sign, word* integer);
void Addition(bigint_st* A, bigint_st* B, bigint_st* C);
void Subtraction(bigint_st* A, bigint_st* B, bigint_st* C);
void Swap_Two_Integers(bigint_st* NUMBER, bigint_st* NUMBER2);
int Compare(bigint_st* A, bigint_st* B);

void Addition_in_Fp(bigint_st* C, bigint_st* P, bigint_st* out);
void Subtraction_in_Fp(bigint_st* C, bigint_st* P, bigint_st* out);
void File_Read(FILE* fp, FILE* fp2, bigint_st* A, bigint_st* B);
void File_Write(FILE* fp, bigint_st* NUMBER);
void Multiplication_PS_64bit(bigint_st* A, bigint_st* B, bigint_st_mul* C);
void Show_integer_Mul(bigint_st_mul* P);
void Set_integer_Mul(bigint_st_mul* P,int sign,  word* integer);
void Multiplication_OS_64bit(bigint_st* A, bigint_st* B, bigint_st* C);
void PS_64_Multiplication(bigint_st* bi_X, bigint_st* bi_Y, bigint_st* bi_Z);
void OS_64_Multiplication(bigint_st* bi_X, bigint_st* bi_Y, bigint_st* bi_Z);
void Multiplication_32bit_Core(word A, word B, word* C);
void Multiplication_Core_64bit(word A, word B, unsigned long long* C);

void word_copy(word* dst, word* src);
void Multiplication_OS_32bit(bigint_st* A, bigint_st* B, bigint_st* C);
void Multiplication_PS_32bit(bigint_st* A, bigint_st* B, bigint_st* C);
void File_Write_Mul(FILE* fp, bigint_st_mul* NUMBER);
void Squaring_64bit(bigint_st* A, bigint_st_mul* C);

void Fast_Reduction(bigint_st_mul* C, bigint_st* P, bigint_st* out);
void integer_copy(bigint_st* bi_dst, bigint_st* bi_src);

void Fermat_Based_Inversion(bigint_st* Z, bigint_st* P, bigint_st* out);

void Elliptic_Curve_Add(ECC_st* P, ECC_st* Q, ECC_st* R, bigint_st* Prime);
void Elliptic_Curve_DBL(ECC_st* P,  ECC_st* R, bigint_st* Prime);

void Set_ECC(ECC_st* P, int sign,  word* integer, word* integer2);
void Show_ECC(ECC_st* P);
void ECC_copy(ECC_st* dst, ECC_st* src);
void Left_to_Right(ECC_st* Base, ECC_st* EC_R, bigint_st* Scalar, bigint_st* Prime);
void File_Write_ECC(FILE* fp, ECC_st* NUMBER);
void Right_to_Left(ECC_st* EC_Base, ECC_st* EC_R, bigint_st* Scalar, bigint_st* Prime);

void Binary_Inversion(bigint_st* A, bigint_st* C, bigint_st* P);
int Word_Compare(bigint_st* A, bigint_st* B);
void Divide_2(bigint_st* NUMBER);

int IsZero(bigint_st* A);

void Jacobian_DBL(ECC_Jacobian_st* P, ECC_Jacobian_st* R, bigint_st* Prime);

void Affine_To_Jacobian(ECC_st* P, ECC_Jacobian_st* R, bigint_st* Prime);

void Jacobian_To_Affine(ECC_Jacobian_st* P, ECC_st* R, bigint_st* Prime);

void Jacobian_ADD(ECC_Jacobian_st* P, ECC_st* Q, ECC_Jacobian_st* R, bigint_st* Prime);

void Jacobian_Copy(ECC_Jacobian_st* dst, ECC_Jacobian_st* src);
void Jacobian_Left_to_Right(ECC_Jacobian_st* EC_Base, ECC_Jacobian_st* EC_R, bigint_st* Scalar, bigint_st* Prime);

void File_Write_Jacobian(FILE* fp, ECC_Jacobian_st* NUMBER);

void Show_ECC_Jacobian(ECC_Jacobian_st* P);


word modular_16(word B);
void NAFw_Fuction(bigint_st* A, char* B,bigint_st* Prime);

void Compute_Table(ECC_st* EC_Base, ECC_st A[],bigint_st* Prime);
void NAF_Left_To_Right_Affine(ECC_st* EC_Base, ECC_st* EC_R, bigint_st* Scalar, bigint_st* Prime); // NAF Affine 버전입니다.

void NAF_Left_To_Right_Jacobian(ECC_Jacobian_st* EC_Base, ECC_Jacobian_st* EC_R, bigint_st* Scalar, bigint_st* Prime);
void Compute_Table_Jacobian(ECC_st* EC_Base, ECC_st A[], bigint_st* Prime);

void Comb_Scalar_Multiplication(ECC_Jacobian_st* EC_Base, ECC_Jacobian_st* EC_R, ECC_st Table[], unsigned char T[][32],bigint_st* Scalar, bigint_st* Prime);

void Comb_Table(ECC_st* EC_Base, ECC_st A[],unsigned char Table[][32], bigint_st* Scalar, bigint_st* Prime);



