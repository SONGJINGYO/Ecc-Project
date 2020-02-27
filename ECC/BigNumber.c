#include"ECC.h"


void Set_integer(bigint_st* P, int sign, word* integer) {
	int j = 0;

	for (int j = 0; j < ECC_256; j++)
	{
		(P)->a[j] = integer[j];
	}
	P->carry = 0;
	P->sign = sign;


	return;
}
void Set_integer_Mul(bigint_st_mul* P,int sign,  word* integer) {
	int j = 0;

	for (int j = 0; j < ECC_256+ECC_256; j++)
	{
		(P)->a[j] = integer[j];
	}
	P->carry = 0;
	P->sign = sign;


	return;
}


void Show_integer_Mul(bigint_st_mul* P) {

	int cnt_i = 0;
	
	printf("0x");

	for (cnt_i = ECC_256+ECC_256 - 1; cnt_i >= 0; cnt_i--) // 높은 자리부터 출력합니다.
	{

		printf("%X", P->a[cnt_i]);
	}
	printf("\n");
	return;
}

int Compare(bigint_st* A, bigint_st* B) {
	int cnt_i = 0;
	if (A->carry == 1) {
		return 1;
	}

	else {
		for (cnt_i = ECC_256- 1; cnt_i >= 0; cnt_i--) {

			if (A->a[cnt_i] > B->a[cnt_i]) {
				return 1; // A가 더 큰 경우
			}
			else if (A->a[cnt_i] < B->a[cnt_i]) {
				return 0; // B가 더 큰 경우
			}
		

		}

		return -1; // 같은 경우
	}

}

int Word_Compare(bigint_st* A, bigint_st* B) {
	int cnt_i = 0;

	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {

		if (A->a[cnt_i] > B->a[cnt_i]) {
			return 1; // A가 더 큰 경우
		}
		else if (A->a[cnt_i] < B->a[cnt_i]) {
			return 0; // B가 더 큰 경우
		}

	}
		return -1; // 같은 경우
	

}




void Swap_Two_Integers(bigint_st* NUMBER, bigint_st* NUMBER2) {
	bigint_st temp;
	temp = *NUMBER2;
	*NUMBER2 = *NUMBER;
	*NUMBER = temp;

}
void word_copy(word* dst, word* src) {
	int cnt_i = 0;

	for (cnt_i = 0; cnt_i < ECC_256; cnt_i++) {
		dst[cnt_i] = src[cnt_i];
	}

	return;



}


void integer_copy(bigint_st* bi_dst, bigint_st* bi_src) //? bigint_src를 bigint_dst에 복사하는 함수
{

	for (int cnt_i = 0; cnt_i < ECC_256; cnt_i++)
	{
		(bi_dst)->a[cnt_i] = (bi_src)->a[cnt_i]; // 할당된 워드에 복사할 값을 대입
	}

	(bi_dst)->carry = (bi_src)->carry; //src의 wordlen 맞춰주기
	(bi_dst)->sign = (bi_src)->sign; //src의 wordlen 맞춰주기
}




void ECC_copy(ECC_st* dst, ECC_st* src) //? bigint_src를 bigint_dst에 복사하는 함수
{
	word_copy(((dst)->x.a),( src->x.a));

	(dst)->x.carry = (src)->x.carry;
	(dst)->x.sign = (src)->x.sign;

	word_copy(((dst)->y.a), (src->y.a));

	(dst)->y.carry = (src)->y.carry;
	(dst)->y.sign = (src)->y.sign;


	return;
}

void Show_ECC(ECC_st* P) {

	int cnt_i = 0;
	printf("X=");
	printf("0x");

	for (cnt_i =ECC_256-1; cnt_i >= 0; cnt_i--) // 높은 자리부터 출력합니다.
	{

		printf("%X ", P->x.a[cnt_i]);
	}
	printf("\n");
	printf("Y=");
	printf("0x");

	for (cnt_i = ECC_256-1; cnt_i >= 0; cnt_i--) // 높은 자리부터 출력합니다.
	{

		printf("%X ", P->y.a[cnt_i]);
	}
	printf("\n");
	return;
}


void Show_ECC_Jacobian(ECC_Jacobian_st* P) {

	int cnt_i = 0;
	printf("X=");
	printf("0x");

	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) // 높은 자리부터 출력합니다.
	{

		printf("%X ", P->x.a[cnt_i]);
	}
	printf("\n");
	printf("Y=");
	printf("0x");

	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) // 높은 자리부터 출력합니다.
	{

		printf("%X ", P->y.a[cnt_i]);
	}
	printf("\n");

	printf("Z=");
	printf("0x");
	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) // 높은 자리부터 출력합니다.
	{

		printf("%X ", P->z.a[cnt_i]);
	}
	printf("\n");

	return;
}

void Show_integer(bigint_st* P) {

	int cnt_i = 0;

	printf("0x");

	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) // 높은 자리부터 출력합니다.
	{

		printf("%X ", P->a[cnt_i]);
	}
	printf("\n");
	return;
}

void Set_ECC(ECC_st* P, int sign, word* integer,word* integer2) {
	int j = 0;

	(P)->x.carry = 0;
	(P)->x.sign = sign;

	word_copy((P)->x.a, integer);


	(P)->y.carry = 0;
	(P)->y.sign = sign;

	word_copy((P)->y.a, integer);

	return;
}