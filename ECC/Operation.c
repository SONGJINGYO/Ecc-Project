#include "ECC.h"

void Addition(bigint_st* A, bigint_st* B, bigint_st* C) {
	int cnt_i = 0;
	(C)->a[0] = (A)->a[0] + (B)->a[0];

	for (cnt_i = 1; cnt_i < ECC_256; cnt_i++) {\
		if ((A)->a[cnt_i - 1] > (C)->a[cnt_i - 1]) {
			C->carry = 1;
		}
		else if ((A)->a[cnt_i - 1] < (C)->a[cnt_i - 1]) {
			C->carry = 0;
		}
		
		(C)->a[cnt_i] = (A)->a[cnt_i] + (B)->a[cnt_i] + (C)->carry;
	}
	if ((A)->a[cnt_i - 1] > (C)->a[cnt_i - 1]) {
		C->carry = 1;
	}
	else if ((A)->a[cnt_i - 1] < (C)->a[cnt_i - 1]) {
		C->carry = 0;
	}
	

	return;
}


void Addition_in_Fp(bigint_st* C, bigint_st* P, bigint_st* out) {


	if (((C)->carry == 1) || (Compare(C, P) != 0)) {

		Subtraction(C, P, out);
		C->carry = 0;
		out->carry = 0;

	}
	else {

		C->carry = 0;
		integer_copy(out, C);

	}

	return;
}

void Subtraction_in_Fp(bigint_st* C, bigint_st* P, bigint_st* out) {

	if ((C)->carry == 1) { //음수
	
		Subtraction(P,C, out);
		C->sign = 0;
		C->carry = 0;
		out->sign = 0;
		out->carry = 0;
	}
	else { // 캐리 여기에 추가하고 리덕션에서 지움 이건 되는지 모름

		C->sign = 0;
		C->carry = 0;
		integer_copy(out, C);
	
	}
	return;
}

void Subtraction(bigint_st* A, bigint_st* B, bigint_st* C) {
	int b_prime = 0;
	int b_d_prime = 0;
	int cnt_i = 0;
	int Borrow = 0;
	int T = 0;
	if (Compare(A, B) == 0) {
		Swap_Two_Integers(A, B);
		(C)->sign = 1;
	}
	for (cnt_i = 0; cnt_i < ECC_256; cnt_i++) {
		b_prime = (A)->a[cnt_i] >= Borrow ? 0 : 1; // A의 배열 값과 b의 배열값을 비교해 b_prime값에 대입해줍니다
		T = (A)->a[cnt_i] - Borrow; //clear
		b_d_prime = T >= (B)->a[cnt_i] ? 0 : 1;//b''의 값을 정해줍니다.
		(C)->a[cnt_i] = T - (B)->a[cnt_i]; // 빼기를 결과값 배열에 저장해줍니다.
		Borrow = b_prime + b_d_prime;
	}
	if ((C)->sign == 1) {
		(C)->carry = 1;
		Swap_Two_Integers(A, B);

	}
	else {
		(C)->carry = Borrow;
	}

	return;
}


void Multiplication_32bit_Core(word A, word B, word* C) {// 단일 워드 곱셈 인풋 한워드 아웃풋 두워드

	int carry = 0;
	word A1B1WA0B0[2] = { 0x00, };
	word T0[2] = { 0x00, };
	word T1[2] = { 0x00, };
	word state[2] = { 0x00, };

	A1B1WA0B0[0] = (A & 0x0000ffff) * (B & 0x0000ffff);
	A1B1WA0B0[1] = (A >> 16) * (B >> 16);

	T0[0] = ((A & 0x0000ffff) * ((B >> 16))) << 16;
	T0[1] = ((A & 0x0000ffff) * ((B >> 16))) >> 16;

	state[0] = T0[0] + A1B1WA0B0[0];
	carry = A1B1WA0B0[0] > state[0] ? 1 : 0;
	state[1] = T0[1] + A1B1WA0B0[1] + carry;

	T1[0] = ((A >> 16) * ((B & 0x0000ffff))) << 16;
	T1[1] = ((A >> 16) * ((B & 0x0000ffff))) >> 16;

	C[0] = T1[0] + state[0];
	carry = state[0] > C[0] ? 1 : 0;
	C[1] = T1[1] + state[1] + carry;

	return;
}






void Multiplication_OS_32bit(bigint_st* A, bigint_st* B, bigint_st* C) {

	int cnt_i = 0;
	int cnt_j = 0;
	word UV[2] = { 0x00, };
	word U = 0x00;
	int carry = 0;
	int carry2 = 0;
	int cnt_k = 0;
	for (cnt_i = 0; cnt_i < ECC_256; cnt_i++) {
		U = 0;
		for (cnt_j = 0; cnt_j < ECC_256; cnt_j++) {

			Multiplication_32bit_Core((A)->a[cnt_i], (B)->a[cnt_j], UV); // 결과는 UV에 들어감

			(C)->a[cnt_i + cnt_j] += U; // C[i+j]=C[i+j]+U

			carry = U > (C)->a[cnt_i + cnt_j] ? 1 : 0; //

			UV[0] += (C)->a[cnt_i + cnt_j];

			carry2 = (C)->a[cnt_i + cnt_j] > UV[0] ? 1 : 0;

			UV[1] += carry + carry2;

			U = UV[1];

			(C)->a[cnt_i + cnt_j] = UV[0];

		}
		(C)->a[cnt_i + ECC_256] = UV[1];

	}

	return;
}


void Multiplication_OS_64bit(bigint_st* A, bigint_st* B, bigint_st* C) {

	int cnt_i = 0;
	int cnt_j = 0;
	word U = 0x00;
	int carry = 0;
	int carry2 = 0;
	int cnt_k = 0;
	unsigned long long UV2 = 0x00;
	for (cnt_i = 0; cnt_i < ECC_256; cnt_i++) {
		U = 0;
		for (cnt_j = 0; cnt_j < ECC_256; cnt_j++) {
			UV2 = (unsigned long long)(A)->a[cnt_i] *
				(unsigned long long)(B)->a[cnt_j] + (C)->a[cnt_i + cnt_j] + U;
			U = UV2 >> 32;
			(C)->a[cnt_i + cnt_j] = UV2 & 0x00000000ffffffff;
		}
		(C)->a[cnt_i + ECC_256] = UV2 >> 32;
	}
	return;
}

void OS_64_Multiplication(bigint_st* bi_X, bigint_st* bi_Y, bigint_st* bi_Z)
{
	int cnt_i, cnt_j; //for loop counting variable
	unsigned long long UV = 0x00LL;
	unsigned int U, V = 0x00;
	word result[ECC_256* 2] = { 0x00 };

	for (cnt_i = 0; cnt_i < ECC_256; cnt_i++)
	{
		UV &= 0x00000000ffffffff;
		for (cnt_j = 0; cnt_j < ECC_256; cnt_j++)
		{
			U = UV >> 32;
			UV = result[cnt_i + cnt_j] + ((unsigned long long)bi_X->a[cnt_i] * (unsigned long long)bi_Y->a[cnt_j]) + U;
			V = UV & 0x00000000ffffffff;
			result[cnt_i + cnt_j] = V;
		}
		U = UV >> 32;
		result[cnt_i +ECC_256] = U;
	}
	memcpy(bi_Z->a, result, (ECC_256 + ECC_256) * sizeof(word));
	return;
}
void PS_64_Multiplication(bigint_st* bi_X, bigint_st* bi_Y, bigint_st* bi_Z)
{

	int cnt_i, cnt_j; //for loop counting variable
	word AH = 0x00;
	word AL = 0x00;
	word BH = 0x00;
	word BL = 0x00;
	word tmp = 0x00;
	word temp1[16] = { 0x00 };
	word temp2[16] = { 0x00 };
	word temp3[16] = { 0x00 };
	int carry = 0;

	for (cnt_i = 0; cnt_i < ECC_256; cnt_i++)
	{
		for (cnt_j = 0; cnt_j < ECC_256; cnt_j++)
		{

			AH = (bi_X)->a[cnt_i] >> 16;
			AL = (bi_X)->a[cnt_i] & 0x0000ffff;
			BH = (bi_Y)->a[cnt_j] >> 16;
			BL = (bi_Y)->a[cnt_j] & 0x0000ffff;

			tmp = AL * BL;
			temp1[cnt_i + cnt_j] += tmp;
			carry = temp1[cnt_i + cnt_j] < (tmp) ? 1 : 0;
			temp1[cnt_i + cnt_j + 1] += carry;

			tmp = AH * BH;
			temp1[cnt_i + cnt_j + 1] += tmp;
			carry = temp1[cnt_i + cnt_j + 1] < (tmp) ? 1 : 0;
			temp1[cnt_i + cnt_j + 2] += carry;
		}
	}
	//!
	for (cnt_i = 0; cnt_i < ECC_256; cnt_i++)
	{
		for (cnt_j = 0; cnt_j < ECC_256; cnt_j++)
		{
			AH = (bi_X)->a[cnt_i] >> 16;
			BH = (bi_Y)->a[cnt_j] >> 16;
			AL = (bi_X)->a[cnt_i] & 0x0000ffff;
			BL = (bi_Y)->a[cnt_j] & 0x0000ffff;

			tmp = AL * BH;
			temp2[cnt_i + cnt_j] += tmp;
			carry = temp2[cnt_i + cnt_j] < tmp ? 1 : 0;
			temp2[cnt_i + cnt_j + 1] += carry;

			tmp = AH * BL;
			temp2[cnt_i + cnt_j] += tmp;
			carry = temp2[cnt_i + cnt_j] < tmp ? 1 : 0;
			temp2[cnt_i + cnt_j + 1] += carry;
		}
	}
	for (cnt_i = 0; cnt_i < 2 * ECC_256; cnt_i++)
	{
		if (cnt_i == 0)
		{
			AH = temp2[cnt_i] >> 16;
			temp2[cnt_i] = temp2[cnt_i] << 16;
			continue;
		}
		AL = temp2[cnt_i] >> 16;
		temp2[cnt_i] = temp2[cnt_i] << 16;
		temp2[cnt_i] &= 0xffff0000;
		temp2[cnt_i] ^= AH;
		AH = AL;
	}
	carry = 0;
	for (cnt_i = 0; cnt_i < ECC_256 * 2; cnt_i++) //둘의 WORD_LEN이 같으므로 한번에 계산 가능하다
	{
		temp3[cnt_i] = temp1[cnt_i] + temp2[cnt_i] + carry; // 단순 덧셈. modulo는 자동적으로 작동
		carry = temp3[cnt_i] < temp1[cnt_i] ? 1 : 0;
		(bi_Z)->a[cnt_i] = temp3[cnt_i];
	}
}

void Multiplication_Core_64bit(word A, word B, unsigned long long* C) {
	unsigned long long state = 0x00;
	unsigned long long state2 = 0x00;

	state ^= (A >> 16) * (B >> 16);
	state = state << 32;
	state ^= (A & 0x0000ffff) * (B & 0x0000ffff);


	state2 ^= ((A & 0x0000ffff) * ((B >> 16))) >> 16;
	state2 = state2 << 32;
	state2 ^= ((A & 0x0000ffff) * ((B >> 16))) << 16;

	state += state2;


	state2 = ((A >> 16) * ((B & 0x0000ffff))) >> 16;
	state2 = state2 << 32;
	state2 ^= ((A >> 16) * ((B & 0x0000ffff))) << 16;
	state += state2;


	*C = state;
	return;

}


void Multiplication_PS_64bit(bigint_st* A, bigint_st* B, bigint_st_mul* C) {

	word R0 = 0x00;
	word R1 = 0x00;
	word R2 = 0x00;
	word V = 0x00;
	word U = 0x00;
	int cnt_i = 0;
	int carry = 0x00;
	int carry2 = 0x00;
	int cnt_j = 0;
	int cnt_k = 0;
	int index = 1;
	unsigned long long UV = 0;
	word temp[16] = { 0x00, };

	for (cnt_k = 0; cnt_k < 2 * ECC_256 - 1; cnt_k++) {
		
		for (cnt_i = cnt_k, cnt_j = 0; cnt_i + cnt_j == cnt_k; cnt_i--, cnt_j++) {
			if (cnt_i >= 8) {
				cnt_i -= index;
				cnt_j += index;
				index++;
			}
			UV = (unsigned long long)(A)->a[cnt_i] * (unsigned long long)(B)->a[cnt_j];

			V = UV & 0x00000000ffffffff;
			U = UV >> 32;

			R0 += V;

			if (V > R0) {
				carry = 1;
			}
			else if (V < R0) {
				carry = 0;
			}
			
			R1 = R1 + U + carry;


			if (U > R1) {
				carry2 = 1;
			}
			else if (U < R1) {
				carry2 = 0;
			}
	
			R2 += carry2;
	
			if (cnt_i == 0 || cnt_j >= 7) {
				break;
	
			}
		
		}
	
		(C)->a[cnt_k] = R0;
	
		R0 = R1;
	
		R1 = R2;
		
		R2 = 0;
	
	}

	(C)->a[(2 * ECC_256) - 1] = R0;

	return;
}






void Multiplication_PS_32bit(bigint_st* A, bigint_st* B, bigint_st* C) {
	word UV[2] = { 0x00, };
	word R0 = 0x00;
	word R1 = 0x00;
	word R2 = 0x00;
	word V = 0x00;
	word U = 0x00;
	int cnt_i = 0;
	int carry = 0x00;
	int carry2 = 0x00;
	int cnt_j = 0;
	int cnt_k = 0;
	int index = 1;
	int cnt_l = 0;
	for (cnt_k = 0; cnt_k < 2 * ECC_256 - 1; cnt_k++) {
		for (cnt_i = cnt_k, cnt_j = 0; cnt_i + cnt_j == cnt_k; cnt_i--, cnt_j++) {
			if (cnt_i >= 8) {
				cnt_i -= index;
				cnt_j += index;
				index++;
			}
			Multiplication_32bit_Core((A)->a[cnt_i], (B)->a[cnt_j], UV); //부하가 큼...
			V = UV[0];
			U = UV[1];
			R0 += V;
			carry = V > R0 ? 1 : 0;
			R1 = R1 + U + carry;
			carry2 = U > R1 ? 1 : 0;
			R2 += carry2;
			if (cnt_i == 0 || cnt_j >= 7) {
				break;
			}

		}
		(C)->a[cnt_k] = R0;
		R0 = R1;
		R1 = R2;
		R2 = 0;
	}
	(C)->a[15] = R0;
	return;
}


void Fast_Reduction(bigint_st_mul* C, bigint_st* P, bigint_st* out) {
	word S1[ECC_256] = { 0x00, };
	word S2[ECC_256] = { 0x00, };
	word S3[ECC_256] = { 0x00, };
	word S4[ECC_256] = { 0x00, };
	word S5[ECC_256] = { 0x00, };
	word S6[ECC_256] = { 0x00, };
	word S7[ECC_256] = { 0x00, };
	word S8[ECC_256] = { 0x00, };
	word S9[ECC_256] = { 0x00, };
	word temp_p[ECC_256] = { 0x00, };
	int cnt_i = 0;
	bigint_st temp = { {0x00},0x00 };
	bigint_st BN_S1 = { {0x00},0x00 };
	bigint_st BN_S2 = { {0x00},0x00 };
	bigint_st BN_S3 = { {0x00},0x00 };



	for (cnt_i = 0; cnt_i < ECC_256; cnt_i++) {
		S1[cnt_i] = (C)->a[cnt_i];
	}
	S2[7] = (C)->a[15]; S2[6] = (C)->a[14]; S2[5] = (C)->a[13]; S2[4] = (C)->a[12]; S2[3] = (C)->a[11];
	S3[6] = (C)->a[15]; S3[5] = (C)->a[14]; S3[4] = (C)->a[13]; S3[3] = (C)->a[12];
	S4[7] = (C)->a[15]; S4[6] = (C)->a[14]; S4[2] = (C)->a[10]; S4[1] = (C)->a[9]; S4[0] = (C)->a[8];
	S5[7] = (C)->a[8]; S5[6] = (C)->a[13]; S5[5] = (C)->a[15]; S5[4] = (C)->a[14]; S5[3] = (C)->a[13]; S5[2] = (C)->a[11]; S5[1] = (C)->a[10]; S5[0] = (C)->a[9];
	S6[7] = (C)->a[10]; S6[6] = (C)->a[8]; S6[2] = (C)->a[13]; S6[1] = (C)->a[12]; S6[0] = (C)->a[11];
	S7[7] = (C)->a[11]; S7[6] = (C)->a[9]; S7[3] = (C)->a[15]; S7[2] = (C)->a[14]; S7[1] = (C)->a[13]; S7[0] = (C)->a[12];
	S8[7] = (C)->a[12]; S8[5] = (C)->a[10]; S8[4] = (C)->a[9]; S8[3] = (C)->a[8]; S8[2] = (C)->a[15]; S8[1] = (C)->a[14]; S8[0] = (C)->a[13];
	S9[7] = (C)->a[13]; S9[5] = (C)->a[11]; S9[4] = (C)->a[10]; S9[3] = (C)->a[9]; S9[1] = (C)->a[15]; S9[0] = (C)->a[14];

	Set_integer(&temp, 0, temp_p);
	Set_integer(&BN_S1, 0, S1);
	Set_integer(&BN_S2, 0, S2);
	Set_integer(&BN_S3, 0, S3);





	Addition(&BN_S1, &BN_S2, &temp);
	Addition_in_Fp(&temp, P, out);


	Addition(out, &BN_S2, &temp);
	Addition_in_Fp(&temp, P, out);


	Addition(out, &BN_S3, &temp);
	Addition_in_Fp(&temp, P, out);


	Addition(out, &BN_S3, &temp);
	Addition_in_Fp(&temp, P, out);



	word_copy(BN_S1.a, S4);
	Addition(out, &BN_S1, &temp);
	Addition_in_Fp(&temp, P, out);


	word_copy(BN_S1.a, S5);
	Addition(out, &BN_S1, &temp);
	Addition_in_Fp(&temp, P, out);

	out->carry = 0;


	word_copy(BN_S1.a, S6);
	Subtraction(out, &BN_S1, &temp);
	Subtraction_in_Fp(&temp, P, out);

	word_copy(BN_S1.a, S7);
	
	Subtraction(out, &BN_S1, &temp);

	Subtraction_in_Fp(&temp, P, out);

	word_copy(BN_S1.a, S8);
	Subtraction(out, &BN_S1, &temp);
	Subtraction_in_Fp(&temp, P, out);




	word_copy(BN_S1.a, S9);
	Subtraction(out, &BN_S1, &temp);
	Subtraction_in_Fp(&temp, P, out);




	return;


}


void Squaring_64bit(bigint_st* A, bigint_st_mul* C) { //A의 워드랜은 8 C의 워드랜은 16

	int cnt_i = 0;
	int cnt_j = 0;
	int cnt_k = 0;
	int index = 1;
	word R1 = 0x00;
	word R0 = 0x00;
	word R2 = 0x00;
	int carry = 0x00;
	word U = 0x00;
	word V = 0x00;
	int temp = 0x00;
	unsigned long long UV = 0;
	for (cnt_k = 0; cnt_k < 2 * ECC_256 - 1; cnt_k++) {
		for (cnt_i = cnt_k, cnt_j = 0; cnt_i + cnt_j == cnt_k; cnt_i--, cnt_j++) {
			if (cnt_i >= 8) {
				cnt_i -= index;
				cnt_j += index;
				index++;
			}

			UV = (unsigned long long)(A)->a[cnt_i] * (unsigned long long)(A)->a[cnt_j];

			V = UV & 0x00000000ffffffff;
			U = UV >> 32;

			if (cnt_i < cnt_j) {
				if ((U & 0x80000000) == 1) {
					UV = 2 * UV;
					carry = 1;
					R2 += carry;
				}
				else {
					UV = 2 * UV;
					carry = 0;
				}
			}

			R0 += V;
			carry = V > R0 ? 1 : 0;
			R1 = R1 + U + carry;
			carry = U > R1 ? 1 : 0;

			R2 += carry;


			if (cnt_i == 0 || cnt_j >= 7) {
				break;
			}


		}(C)->a[cnt_k] = R0;
		R0 = R1;
		R1 = R2;
		R2 = 0;



	}(C)->a[2 * ECC_256 - 1] = R0;
	return;
}

void Fermat_Based_Inversion(bigint_st* Z, bigint_st* P, bigint_st* out) {
	bigint_st Copy_Integer = { {0x00},0x00 };
	bigint_st_mul Mul_temp = { {0x00},0x00 };
	bigint_st T0 = { {0x00},0x00 };
	bigint_st T1 = { {0x00},0x00 };
	int cnt_i = 0;
	int LtR = 0;
	int temp = 0;
	word Mul_temp_p[ECC_256 + ECC_256] = { 0x00, };
	word temp_p[ECC_256] = { 0x00, };
	Set_integer(&T0, 0,  temp_p);
	Set_integer(&T1, 0,  temp_p);
	Set_integer_Mul(&Mul_temp, 0,  Mul_temp_p);
	Set_integer(&Copy_Integer, 0, temp_p);


	//z^3
	Squaring_64bit(Z, &Mul_temp); // z^2을 저장

	Fast_Reduction(&Mul_temp, P, out); // z^2 Reduction out에 저장

	Multiplication_PS_64bit(out, Z, &Mul_temp); // z^2*z 곱하기

	Fast_Reduction(&Mul_temp, P, out); // z^3 리덕션 결과 out에 저장

	integer_copy(&Copy_Integer, out); // z^3을 저장

	//z^15
	Squaring_64bit(out, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	Squaring_64bit(out, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out); // z^12승 리덕션까지 적용

	Multiplication_PS_64bit(out, &Copy_Integer, &Mul_temp); //z^15 을 Mul_temp에 저장

	Fast_Reduction(&Mul_temp, P, out); //리덕션 적용

	//t0 계산이 필요함
	Squaring_64bit(out, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	Squaring_64bit(out, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	Multiplication_PS_64bit(out, &Copy_Integer, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	integer_copy(&T0, out);

	//t1 계산이 필요함

	for (cnt_i = 0; cnt_i < 6; cnt_i++) {
		Squaring_64bit(out, &Mul_temp);

		Fast_Reduction(&Mul_temp, P, out);

	}

	Multiplication_PS_64bit(out, &T0, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	integer_copy(&T1, out); // T1을 저장 

	//t2 계산이 필요함

	for (cnt_i = 0; cnt_i < 12; cnt_i++) {
		Squaring_64bit(out, &Mul_temp);

		Fast_Reduction(&Mul_temp, P, out);

	}

	Multiplication_PS_64bit(out, &T1, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);


	for (cnt_i = 0; cnt_i < 6; cnt_i++) {
		Squaring_64bit(out, &Mul_temp);

		Fast_Reduction(&Mul_temp, P, out);

	}

	Multiplication_PS_64bit(out, &T0, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	integer_copy(&T0, out); // T2를 저장

	//t3 계산 필요

	Squaring_64bit(out, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	Squaring_64bit(out, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);


	Multiplication_PS_64bit(out, &Copy_Integer, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	integer_copy(&T1, out);//t3을 저장시킴 t1에

	//t4 계산 필요 여기 지우면 됨...
	for (cnt_i = 0; cnt_i < 32; cnt_i++) {
		Squaring_64bit(out, &Mul_temp);

		Fast_Reduction(&Mul_temp, P, out);

	}

	Multiplication_PS_64bit(out, Z, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	for (cnt_i = 0; cnt_i < 96; cnt_i++) {
		Squaring_64bit(out, &Mul_temp);
		Fast_Reduction(&Mul_temp, P, out);

	}


	//t5 계산

	for (cnt_i = 0; cnt_i < 32; cnt_i++) {
		Squaring_64bit(out, &Mul_temp);

		Fast_Reduction(&Mul_temp, P, out);

	}

	Multiplication_PS_64bit(out, &T1, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	for (cnt_i = 0; cnt_i < 32; cnt_i++) {
		Squaring_64bit(out, &Mul_temp);

		Fast_Reduction(&Mul_temp, P, out);

	}
	Multiplication_PS_64bit(out, &T1, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	//t 계산
	for (cnt_i = 0; cnt_i < 30; cnt_i++) {
		Squaring_64bit(out, &Mul_temp);

		Fast_Reduction(&Mul_temp, P, out);

	}
	Multiplication_PS_64bit(out, &T0, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);


	Squaring_64bit(out, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);

	Squaring_64bit(out, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);
	Multiplication_PS_64bit(out, Z, &Mul_temp);

	Fast_Reduction(&Mul_temp, P, out);


	return;
}



void Elliptic_Curve_Add(ECC_st* P, ECC_st* Q, ECC_st* R, bigint_st* Prime) { // P,Q 더렵히지지만 결과값은 나옴 즉 P,Q
	bigint_st temp = { {0x00},0x00 };
	bigint_st temp2 = { {0x00},0x00 };
	bigint_st temp3 = { {0x00},0x00 };

	bigint_st_mul Mul_temp = { {0x00},0x00 };
	word temp_p[ECC_256] = { 0x00, };
	word Mul_temp_p[ECC_256 + ECC_256] = { 0x00, };

	Set_integer(&temp, 0,  temp_p);
	Set_integer(&temp2, 0,  temp_p);
	Set_integer(&temp3, 0,  temp_p);

	Set_integer_Mul(&Mul_temp, 0,Mul_temp_p);

	Subtraction(&((Q)->y), &((P)->y), &temp);
	Subtraction_in_Fp(&temp, Prime, &temp3);   // y2-y1의 결과값을 temp3에 저장


	Subtraction(&((Q)->x), &((P)->x), &temp);
	Subtraction_in_Fp(&temp, Prime, &temp2); // x2-x1 값을 temp2에 저장


	Fermat_Based_Inversion(&temp2, Prime, &temp); // x2-x1 값의 역원을 temp에 저장


	Multiplication_PS_64bit(&temp, &temp3, &Mul_temp);
	Fast_Reduction(&Mul_temp, Prime, &temp); // 기울기 값을 저장


	Squaring_64bit(&temp, &Mul_temp);
	Fast_Reduction(&Mul_temp, Prime, &temp2); // 제곱값의 리덕션 값을 y2에 저장


	Subtraction(&temp2, &((P)->x), &temp3);
	Subtraction_in_Fp(&temp3, Prime, &temp2);  // -x1의 값을 리덕션 적용해서 x1에 저장


	Subtraction(&temp2, &((Q)->x), &temp3);
	Subtraction_in_Fp(&temp3, Prime, &((R)->x)); // 최종 결과값을 R의  x 좌표에 저장

	integer_copy(&temp2, &((R)->x));

	Subtraction(&((P)->x), &temp2, &temp3);
	Subtraction_in_Fp(&temp3, Prime, &temp2); //

	Multiplication_PS_64bit(&temp2, &temp, &Mul_temp);
	Fast_Reduction(&Mul_temp, Prime, &temp); // 기울기 값을 저장

	Subtraction(&temp, &((P)->y), &temp2);
	Subtraction_in_Fp(&temp2, Prime, &((R)->y));

	// 값이 섞였는지 확인이 필요함.
	return;
}




void Elliptic_Curve_DBL(ECC_st* P, ECC_st* R, bigint_st* Prime) { // P,Q 더렵히지지만 결과값은 나옴 즉 P,Q
	bigint_st temp = { {0x00},0x00 };
	bigint_st temp2 = { {0x00},0x00 };
	bigint_st temp3 = { {0x00},0x00 };
	bigint_st temp4 = { {0x00},0x00 };
	bigint_st NUMBER_3 = { {0x00},0x00 };

	bigint_st_mul Mul_temp = { {0x00},0x00 };
	word temp_p[ECC_256] = { 0x00, };
	word temp_p3[ECC_256] = { 3,0x00,0x00,0x00,0x00,0x00,0x00,0x00 };

	word Mul_temp_p[ECC_256 + ECC_256] = { 0x00, };

	Set_integer_Mul(&Mul_temp, 0,  Mul_temp_p);
	Set_integer(&temp, 0, temp_p);
	Set_integer(&NUMBER_3, 0,temp_p3); // 숫자 3
	Set_integer(&temp2, 0, temp_p);
	Set_integer(&temp3, 0,  temp_p);
	Set_integer(&temp4, 0,  temp_p);

	Squaring_64bit(&((P)->x), &Mul_temp);
	Fast_Reduction(&Mul_temp, Prime, &temp); // x1 제곱


	
	Multiplication_PS_64bit(&temp, &NUMBER_3, &Mul_temp);
	Fast_Reduction(&Mul_temp, Prime, &temp); //3을 곱함


	Subtraction(Prime, &NUMBER_3, &temp2);
	Subtraction_in_Fp(&temp2, Prime, &temp4); // 윗줄 끝
	
	
	Addition(&temp4, &temp, &temp2);
	Addition_in_Fp(&temp2, Prime, &temp);
	
	integer_copy(&temp2, &((P)->y));
	Addition(&((P)->y), &temp2, &temp3); //
	Addition_in_Fp(&temp3, Prime, &temp2); //2 y1
	//Okay
	
	Fermat_Based_Inversion(&temp2, Prime, &temp3); /// No touch temp3


	Multiplication_PS_64bit(&temp3, &temp, &Mul_temp);
	Fast_Reduction(&Mul_temp, Prime, &temp3);
	
	//Okay
	Squaring_64bit(&temp3, &Mul_temp);
	Fast_Reduction(&Mul_temp, Prime, &temp); //
	
	//
	integer_copy(&temp2, &(P->x));
	Addition(&(P->x), &temp2, &temp4);
	Addition_in_Fp(&temp4, Prime, &temp2);
	
	Subtraction(&temp, &temp2, &temp4);
	Subtraction_in_Fp(&temp4, Prime, &(R->x));
	


	integer_copy(&temp4, &(R->x));
	Subtraction(&(P->x), &temp4, &temp2);
	Subtraction_in_Fp(&temp2, Prime, &temp);

	Multiplication_PS_64bit(&temp3, &temp, &Mul_temp);
	Fast_Reduction(&Mul_temp, Prime, &temp3);

	Subtraction(&temp3, &(P->y), &temp);
	Subtraction_in_Fp(&temp, Prime, &(R->y));



	// 값이 섞였는지 확인이 필요함.
	return;
}

void Left_to_Right(ECC_st* EC_Base, ECC_st* EC_R, bigint_st* Scalar, bigint_st* Prime) {

	int cnt_i = 0;
	int cnt_j = 0;
	int temp = 0;
	
	ECC_st EC_P= { {0x00},0x00 };
	word P[ECC_256] = { 0x00, };
	Set_ECC(&EC_P, 0, P, P);
	int first_bit = 0;

	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		for (cnt_j = 31; cnt_j >= 0; cnt_j--) {
		
		
			temp = ((Scalar->a[cnt_i]) >> cnt_j) & 0x1;
			
if (temp == 1) {
	first_bit += 1;
	if (first_bit == 1) {
		ECC_copy(&EC_P, EC_Base);
	}
	else if (first_bit > 1) {
		Elliptic_Curve_DBL(&EC_P, EC_R, Prime);
		Elliptic_Curve_Add(EC_R, EC_Base, &EC_P, Prime);
		
	}

}
else if (temp == 0) {
	if (first_bit >= 1) {
		Elliptic_Curve_DBL(&EC_P, EC_R, Prime);
		ECC_copy(&EC_P, EC_R);
	}

}


		}
		ECC_copy(EC_R, &EC_P);
	}
	return;
}
void Right_to_Left(ECC_st* EC_Base, ECC_st* EC_R, bigint_st* Scalar, bigint_st* Prime) {
	int cnt_i = 0;
	int cnt_j = 0;
	int temp = 0;
	ECC_st EC_P = { {0x00},0x00 };
	ECC_st EC_Q = { {0x00},0x00 };
	word P[ECC_256] = { 0x00, };
	Set_ECC(&EC_P, 0, P, P); /// 베이스 스타트 0
	Set_ECC(&EC_Q, 0, P, P);/// 스타트 P
	int first_bit = 0;

	ECC_copy(&EC_Q, EC_Base); // 스타트 P
	for (cnt_i = 0; cnt_i < ECC_256; cnt_i++) {
		for (cnt_j = 0; cnt_j < 32; cnt_j++) {

			temp = ((Scalar->a[cnt_i]) >> cnt_j) & 0x1;

			if (temp == 1) {
				first_bit += 1;
				if (first_bit == 1) {
					ECC_copy(&EC_P, &EC_Q);
					Elliptic_Curve_DBL(&EC_Q, EC_R, Prime);
					ECC_copy(&EC_Q, EC_R);
				}
				else {
					Elliptic_Curve_Add(&EC_P, &EC_Q, EC_R, Prime);
					ECC_copy(&EC_P, EC_R);

					Elliptic_Curve_DBL(&EC_Q, EC_R, Prime);
					ECC_copy(&EC_Q, EC_R);
				}
			}


			else {
				Elliptic_Curve_DBL(&EC_Q, EC_R, Prime);
				ECC_copy(&EC_Q, EC_R);
			}



		}
		ECC_copy(EC_R, &EC_P);
	}
	return;

}
void Divide_2(bigint_st* NUMBER) {
	int cnt_i = 0;
	int temp = 0;
	for (cnt_i = 0; cnt_i < ECC_256; cnt_i++) {
		(NUMBER)->a[cnt_i] = (NUMBER)->a[temp + 1] << 31 | (NUMBER)->a[temp] >> 1;
		temp += 1;
	}
	return;
}



void Binary_Inversion(bigint_st* A, bigint_st* C, bigint_st* P) {
	bigint_st temp = { {0x00},0x00 };
	bigint_st temp2 = { {0x00},0x00 };
	bigint_st temp3= { {0x00},0x00 };
	bigint_st Number_1 = { {0x00},0x00 };
	bigint_st X1 = { {0x00},0x00 };
	bigint_st X2 = { {0x00},0x00 };
	word temp_p[ECC_256] = { 1,0x00,0x00,0x00,0x00,0x00,0x00,0x00 };
	word X2_temp[ECC_256] = { 0x00, };
	Set_integer(&Number_1, 0, temp_p);
	Set_integer(&X2, 0, X2_temp);
	Set_integer(&temp3, 0, X2_temp);


	int cnt_i = 0;

	int cnt_j = 0;
	int cnt_k = 0;
	integer_copy(&temp, A); // U


	integer_copy(&temp2, P); //V

	integer_copy(&X1, &Number_1);


	while ((Word_Compare(&temp, &Number_1) != -1) && (Word_Compare(&temp2, &Number_1) != -1)) {
	
		cnt_i++;

		while ((temp.a[0] & 1) == 0) { // 짝수 인 경우
		
			cnt_j++;
			Divide_2(&temp);

			if ((X1.a[0] & 1) == 0){
				Divide_2(&X1);

			}
			else{
	
				Addition(&X1, P, &temp3);
				if (temp3.carry == 1) {
					Divide_2(&temp3);
					temp3.a[ECC_256 - 1] ^= 0x80000000;
				}
				else {
					Divide_2(&temp3);
				}
				temp3.carry = 0;
		
				Addition_in_Fp(&temp3, P, &X1);
			


	
			}
			
		}

		while ((temp2.a[0] & 1) == 0) { // 짝수 인 경우
			Divide_2(&temp2);
	
			cnt_k++;
			if ((X2.a[0] & 1) == 0) {
				Divide_2(&X2);
			}
			else {

				Addition(&X2, P, &temp3);
				if (temp3.carry == 1) {
					Divide_2(&temp3);
					temp3.a[ECC_256 - 1] ^= 0x80000000;
				}
				else {
					Divide_2(&temp3);
				}

				
				temp3.carry = 0;
				Addition_in_Fp(&temp3, P, &X2);
				

			}
		
		}
			if (Word_Compare(&temp, &temp2) != 0) {
			
				Subtraction(&temp, &temp2, &temp3);

		
				Subtraction_in_Fp(&temp3, P, &temp);

				Subtraction(&X1, &X2, &temp3);
				Subtraction_in_Fp(&temp3, P, &X1);
				

			}
			else {
			

				Subtraction(&temp2, &temp, &temp3);

				Subtraction_in_Fp(&temp3, P, &temp2);		
		
				Subtraction(&X2, &X1, &temp3);
				Subtraction_in_Fp(&temp3, P, &X2);
		
		
			}
	}

	if (Word_Compare(&Number_1, &temp) == -1) {
		
		integer_copy(C, &X1);
	}
	else {
		integer_copy(C, &X2);
	}
	return;

}
void Affine_To_Jacobian(ECC_st* P, ECC_Jacobian_st* R,bigint_st* Prime) {
	bigint_st temp= { {0x00},0x00 };
	bigint_st temp2 = { {0x00},0x00 };
	bigint_st_mul temp3 = { {0x00},0x00 };
	word temp_p[ECC_256] = { 0x00, };
	word temp_p2[ECC_256+ECC_256] = { 0x00, };
	Set_integer(&temp, 0, temp_p);
	Set_integer(&temp2, 0, temp_p);
	Set_integer_Mul(&temp3, 0, temp_p2);
	
	integer_copy(&(R->z), &temp);
	R->z.a[0] = 1;
	R->x = P->x;
	R->y = P->y;
	

	


	return;

}


void Jacobian_To_Affine(ECC_Jacobian_st* P, ECC_st* R, bigint_st* Prime) {

	bigint_st temp = { {0x00},0x00 };

	bigint_st_mul temp2 = { {0x00},0x00 };
	bigint_st temp3 = { {0x00},0x00 };
	word temp_p[ECC_256] = { 0x00, };
	word temp_p2[ECC_256+ECC_256] = { 0x00, };
	Set_integer(&temp, 0, temp_p);
	Set_integer_Mul(&temp2, 0, temp_p2);
	Set_integer(&temp3, 0, temp_p);
	
	integer_copy(&temp3, &(P->z));

	Squaring_64bit(&(P->z), &temp2);
	Fast_Reduction(&temp2, Prime, &(P->z)); // Z^2

	Multiplication_PS_64bit(&(P->z), &temp3, &temp2);
	Fast_Reduction(&temp2, Prime, &(P->z));

	Binary_Inversion(&(P->z), &temp, Prime);

	Multiplication_PS_64bit(&temp, &(P->y), &temp2);
	Fast_Reduction(&temp2, Prime, &(R->y));

	Multiplication_PS_64bit(&temp, &temp3, &temp2);
	Fast_Reduction(&temp2, Prime, &(P->z));

	
	Multiplication_PS_64bit(&(P->z), &(P->x), &temp2);
	Fast_Reduction(&temp2, Prime, &(R->x));





	return;
}



void Jacobian_DBL(ECC_Jacobian_st* P,ECC_Jacobian_st* R,bigint_st* Prime) {
	bigint_st temp= { {0x00},0x00 };
	bigint_st temp2 = { {0x00},0x00 };
	bigint_st T1 = { {0x00},0x00 };
	bigint_st T2 = { {0x00},0x00 };
	bigint_st T3 = { {0x00},0x00 };

	bigint_st_mul mul_temp = { {0x00},0x00 };
	word temp_p[ECC_256] = { 0x00, };
	word temp_mul[ECC_256+ECC_256] = { 0x00, };

	Set_integer(&temp, 0, temp_p);
	Set_integer(&temp2, 0, temp_p);
	Set_integer(&T2, 0, temp_p);
	Set_integer(&T3, 0, temp_p);
	Set_integer_Mul(&mul_temp, 0, temp_mul);
	Set_integer(&T1, 0, temp_p);

	//1
	Squaring_64bit(&(P->z), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T1);

	//2
	Subtraction(&(P->x), &T1, &temp);
	Subtraction_in_Fp(&temp, Prime, &T2);

	//3
	Addition(&(P->x), &T1, &temp);
	Addition_in_Fp(&temp, Prime, &T1);

	//4
	Multiplication_PS_64bit(&T2, &T1, &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T2);

	//5
	integer_copy(&temp2, &T2);

	Addition(&temp2, &T2, &temp);
	Addition_in_Fp(&temp, Prime, &T2);

	Addition(&temp2, &T2, &temp);
	Addition_in_Fp(&temp, Prime, &T2);

	//6
	integer_copy(&temp, &(P->y));

	
	Addition(&(P->y), &temp, &temp2);
	Addition_in_Fp(&temp2, Prime, &(R->y));

	//7
	Multiplication_PS_64bit(&(R->y), &(P->z), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &(R->z));

	//8
	Squaring_64bit(&(R->y), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &(R->y));

	//9
	Multiplication_PS_64bit(&(R->y), &(P->x), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T3);

	//
	Squaring_64bit(&(R->y), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &(R->y));

	if ((R->y.a[0] & 1 )== 0) {
		Divide_2(&(R->y));

	}
	else {
		integer_copy(&temp2, &(R->y));
		Addition(Prime, &temp2, &(R->y));
		if (R->y.carry == 1) {
			Divide_2(&(R->y));
			R->y.a[ECC_256 - 1] ^= 0x80000000;
		}
	}

	Squaring_64bit(&T2, &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &(R->x));

	integer_copy(&temp, &T3);
	Addition(&temp, &T3, &temp2);
	Addition_in_Fp(&temp2, Prime, &T1);

	//15

	Subtraction(&(R->x), &T1, &temp);
	Subtraction_in_Fp(&temp, Prime, &(R->x));
	//16

	Subtraction(&T3,&(R->x), &temp);
	Subtraction_in_Fp(&temp, Prime, &T1);


	Multiplication_PS_64bit(&T1, &T2, &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T1);


	Subtraction(&T1, &(R->y), &temp);
	Subtraction_in_Fp(&temp, Prime, &(R->y));




	return;




}

int IsZero(bigint_st* A) {
	int cnt_i = 0;
	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		if (A->a[cnt_i] != 0) {
			return 1;
		}
	}
	return 0;

}


void Jacobian_ADD(ECC_Jacobian_st* P, ECC_st* Q, ECC_Jacobian_st* R, bigint_st* Prime) {
	bigint_st temp = { {0x00},0x00 };
	bigint_st temp2 = { {0x00},0x00 };
	bigint_st T1 = { {0x00},0x00 };
	bigint_st T2 = { {0x00},0x00 };
	bigint_st T3 = { {0x00},0x00 };
	bigint_st T4 = { {0x00},0x00 };
	ECC_Jacobian_st J_temp= { {0x00},0x00 };

	bigint_st_mul mul_temp = { {0x00},0x00 };
	word temp_p[ECC_256] = { 0x00, };
	word temp_mul[ECC_256 + ECC_256] = { 0x00, };

	Set_integer(&temp, 0, temp_p);
	Set_integer(&temp2, 0, temp_p);
	Set_integer(&T2, 0, temp_p);
	Set_integer(&T3, 0, temp_p);
	Set_integer_Mul(&mul_temp, 0, temp_mul);
	Set_integer(&T1, 0, temp_p);
	Set_integer(&T4, 0, temp_p);
	//1
	Squaring_64bit(&(P->z), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T1);
	//2
	Multiplication_PS_64bit(&T1, &(P->z), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T2);

	//3
	Multiplication_PS_64bit(&T1, &(Q->x), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T1);

	//4
	Multiplication_PS_64bit(&T2, &(Q->y), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T2);
	//5
	Subtraction(&T1, &(P->x), &temp);
	Subtraction_in_Fp(&temp, Prime, &T1);
	//6
	Subtraction(&T2, &(P->y), &temp);
	Subtraction_in_Fp(&temp, Prime, &T2);

	//7
	if ((IsZero(&T1))==0) {
		if (IsZero(&T2) == 0) {
			Affine_To_Jacobian(Q, &J_temp, Prime);
			Jacobian_DBL(&J_temp, R,Prime);
		}
		

	}
	Multiplication_PS_64bit(&(P->z), &T1, &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &(R->z));

	Squaring_64bit(&T1, &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T3);

	Multiplication_PS_64bit(&T3, &T1, &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T4);

	Multiplication_PS_64bit(&T3, &(P->x), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T3);

	integer_copy(&temp, &T3);
	Addition(&temp, &T3, &temp2);
	Addition_in_Fp(&temp2, Prime, &T1);

	Squaring_64bit(&T2, &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &(R->x));

	//16
	Subtraction(&(R->x), &T1, &temp);
	Subtraction_in_Fp(&temp, Prime, &(R->x));

	Subtraction(&(R->x), &T4, &temp);
	Subtraction_in_Fp(&temp, Prime, &(R->x));

	//18

	Subtraction(&T3, &(R->x), &temp);
	Subtraction_in_Fp(&temp, Prime, &T3);

	Multiplication_PS_64bit(&T2, &T3, &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T3);

	Multiplication_PS_64bit(&T4, &(P->y), &mul_temp);
	Fast_Reduction(&mul_temp, Prime, &T4);

	Subtraction(&T3, &T4, &temp);
	Subtraction_in_Fp(&temp, Prime, &(R->y));



	return;


}

void Jacobian_Copy(ECC_Jacobian_st* dst, ECC_Jacobian_st* src) {
	dst->x = src->x;
	dst->y = src->y;
	dst->z = src->z;
	dst->is_infinity = src->is_infinity;
	return;

}
void Jacobian_Left_to_Right(ECC_Jacobian_st* EC_Base, ECC_Jacobian_st* EC_R, bigint_st* Scalar, bigint_st* Prime) {

	int cnt_i = 0;
	int cnt_j = 0;
	int temp = 0;
	ECC_Jacobian_st EC_P = { {0x00},0x00 };
	ECC_st temp2= { {0x00},0x00 };
	word P[ECC_256] = { 0x00, };
	Jacobian_To_Affine(EC_Base, &temp2, Prime);
	
	int first_bit = 0;

	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		for (cnt_j = 31; cnt_j >= 0; cnt_j--) {
		
			temp = ((Scalar->a[cnt_i]) >> cnt_j) & 0x1;

			if (temp == 1) {
				first_bit += 1;
				if (first_bit == 1) {
					Jacobian_Copy(&EC_P, EC_Base);
					
			
				}
				else if (first_bit > 1) {

					Jacobian_DBL(&EC_P, EC_R, Prime);
					
					Jacobian_ADD(EC_R, &temp2, &EC_P, Prime);
				
				}

			}
			else if (temp == 0) {
				if (first_bit >= 1) {
					
					Jacobian_DBL(&EC_P, EC_R, Prime);
					
					Jacobian_Copy(&EC_P, EC_R);

				}

			}


		}
		Jacobian_Copy(EC_R, &EC_P);
	}
	return;
}
word modular_16(word B) {

	word a = 0x0000000f;
	word temp = 0;
	B = a & (B);

	temp ^= B;

	if (temp > 7) {
		temp = temp-16;
	}

	return temp;
}
void NAFw_Fuction(bigint_st* A,char* B,bigint_st* Prime){
	int cnt_i = 0;
	word temp[ECC_256] = { 1,0,0,0,0,0,0,0 };
	word Zero[ECC_256] = { 0x00, };
	bigint_st temp_N= { {0x00},0x00 };
	unsigned char temp3 = 0x00;
	Set_integer(&temp_N, 0, temp);

	bigint_st Sub_temp= { {0x00},0x00 };
	bigint_st Ki= { {0x00},0x00 };
	Set_integer(&Ki, 0, Zero);
	Set_integer(&Sub_temp, 0, Zero);

	while (Word_Compare(A, &temp_N) != 0) {

		if (((A->a[0]) & 1 )== 1) {
	
			B[cnt_i] =modular_16(A->a[0]);
			//Ki.a[0]=(int)B[cnt_i];
			Ki.a[0] = B[cnt_i];
			if (Ki.a[0] > 8) {
				Ki.a[0] = ~(Ki.a[0]) + 1;
				//integer_copy(&Sub_temp, A);
				//Addition(&Sub_temp, &Ki, A);

				Addition(A, &Ki, &Sub_temp);
				Addition_in_Fp(&Sub_temp, Prime, A);

			}
			else {
			//	integer_copy(&Sub_temp, A);
			//	Subtraction(&Sub_temp, &Ki, A);
				Subtraction(A, &Ki, &Sub_temp);
				Subtraction_in_Fp(&Sub_temp, Prime, A);
			}
		
		}
		else {
			B[cnt_i] = 0x00;
		}
		Divide_2(A);
		cnt_i += 1;
	}
	return;

}
void Compute_Table(ECC_st* EC_Base,ECC_st A[],bigint_st* Prime){
	ECC_st Temp= { {0x00},0x00 };
	ECC_st Temp2 = { {0x00},0x00 };
	word temp_ecc[ECC_256] = { 0x00, };
	Set_ECC(&Temp, 0, temp_ecc, temp_ecc);
	Set_ECC(&Temp2, 0, temp_ecc, temp_ecc);


	//1P
	ECC_copy(&A[0], EC_Base);

	//3P
	Elliptic_Curve_DBL(&A[0], &Temp, Prime); // Temp = 2P
	Elliptic_Curve_Add(&Temp, EC_Base, &A[1], Prime); // A[1] = 3P
	//5P
	Elliptic_Curve_Add(&Temp, &A[1], &A[2], Prime); // A[2]=5P
	//7P
	Elliptic_Curve_Add(&A[2], &Temp, &A[3], Prime); // A[2]=5P

	
	
}
void NAF_Left_To_Right_Affine(ECC_st* EC_Base, ECC_st* EC_R, bigint_st* Scalar, bigint_st* Prime) {
	char NAFw[257] = { 0x00, };
	ECC_st Table[4]= { {0x00},0x00 };
	ECC_st EC_Temp={ {0x00},0x00 };
	int cnt_i = 0;
	int cnt_j = 0;
	int temp = 0;
	ECC_st EC_P = { {0x00},0x00 };
	bigint_st BN_Temp = { {0x00},0x00 };
	word P[ECC_256] = { 0x00, };
	
	int first_bit = 0;
	Set_integer(&BN_Temp, 0,P);
	Set_ECC(&EC_P, 0, P, P);

	NAFw_Fuction(Scalar, NAFw, Prime);
	Compute_Table(EC_Base, Table, Prime);
	
	

	for (cnt_i = 256; cnt_i >= 0; cnt_i--) {

		temp = (NAFw[cnt_i]);
		
		if (temp != 0) {
			first_bit += 1;
			if (first_bit == 1) {
				if (temp > 0) {
					ECC_copy(&EC_P, &Table[(temp-1)/2]);
				//	printf("Temp=%d, %d\n", temp, (temp - 1) / 2);
				//	printf("First %d-th temp=%d\n", cnt_i, temp);
				//	Show_ECC(&EC_P);
				//	printf("\n\n");
				}
				else {
					ECC_copy(&EC_P,&Table[(-temp - 1) / 2]);
				//	printf("Temp=%d, %d\n", temp, (-(temp + 1) / 2));
			//		printf("Frist %d-th temp=%d\n", cnt_i,temp);
				//	Show_ECC(&EC_P);
				//	printf("\n\n");
				}

			}
		
			else {
				if (temp > 0) {
					Elliptic_Curve_DBL(&EC_P, EC_R, Prime);
					Elliptic_Curve_Add(EC_R, &Table[(temp-1)/2], &EC_P, Prime);

					//printf("temp>0 %d-th temp=%d, %d\n", cnt_i, temp,(temp-1)/2);
		
				//	Show_ECC(&EC_P);
				//	printf("\n\n");
				}
				else {
					Elliptic_Curve_DBL(&EC_P, EC_R, Prime);
					ECC_copy(&EC_Temp, &Table[(-temp-1)/2]);
					Subtraction(Prime, &(EC_Temp.y), &BN_Temp);
					Subtraction_in_Fp(&BN_Temp, Prime, &(EC_Temp.y));


					Elliptic_Curve_Add(EC_R, &EC_Temp, &EC_P, Prime);
	
				//	printf("temp<0 %d-th temp=%d %d\n", cnt_i, temp,-(temp+1)/2);
	
				//	Show_ECC(&EC_P);
				//	printf("\n\n");
				}

			}
		}
		else {
			if ((first_bit > 0)) {
				Elliptic_Curve_DBL(&EC_P, EC_R, Prime);
				ECC_copy(&EC_P, EC_R);
				//printf("temp=Zero %d-th temp=%d\n", cnt_i, temp);
				//Show_ECC(&EC_P);
			//	printf("\n\n");
			}

		}

		ECC_copy(EC_R, &EC_P);
	}


	return;

}


void Compute_Table_Jacobian(ECC_st* EC_Base, ECC_st A[], bigint_st* Prime) {
	ECC_st Temp = { {0x00},0x00 };
	ECC_st Temp2 = { {0x00},0x00 };
	word temp_ecc[ECC_256] = { 0x00, };
	Set_ECC(&Temp, 0, temp_ecc, temp_ecc);
	Set_ECC(&Temp2, 0, temp_ecc, temp_ecc);

	ECC_Jacobian_st EC_J_Temp= { {0x00},0x00 };

	ECC_Jacobian_st EC_J_Temp3 = { {0x00},0x00 };
	ECC_Jacobian_st EC_J_Temp2= { {0x00},0x00 };

	Affine_To_Jacobian(EC_Base, &EC_J_Temp, Prime); // 자코비안 P 제대로 담김

	
	//1P
	ECC_copy(&A[0], EC_Base);

	//3P
	Jacobian_DBL(&EC_J_Temp, &EC_J_Temp2, Prime); //2P
	Jacobian_Copy(&EC_J_Temp, &EC_J_Temp2);
	Jacobian_To_Affine(&EC_J_Temp2, &Temp, Prime); // Temp에는 아핀으로 2P


	Jacobian_ADD(&EC_J_Temp, EC_Base, &EC_J_Temp2, Prime); //3p
	Jacobian_Copy(&EC_J_Temp, &EC_J_Temp2);
	Jacobian_To_Affine(&EC_J_Temp2, &A[1], Prime); // Temp2에는 2P // Temp에는 3P
	

	//5P
	Jacobian_ADD(&EC_J_Temp, &Temp, &EC_J_Temp2, Prime);
	Jacobian_Copy(&EC_J_Temp, &EC_J_Temp2);
	Jacobian_To_Affine(&EC_J_Temp2, &A[2], Prime); // 5P

	//7P
	Jacobian_ADD(&EC_J_Temp, &Temp, &EC_J_Temp2, Prime);
	Jacobian_Copy(&EC_J_Temp, &EC_J_Temp2);
	Jacobian_To_Affine(&EC_J_Temp2, &A[3], Prime);


	return;

}
void NAF_Left_To_Right_Jacobian(ECC_Jacobian_st* EC_Base, ECC_Jacobian_st* EC_R, bigint_st* Scalar, bigint_st* Prime) {
	char NAFw[257] = { 0x00, };
	ECC_st Table[4] = { {0x00},0x00 };
	ECC_st EC_Temp = { {0x00},0x00 };
	int cnt_i = 0;
	int cnt_j = 0;
	int temp = 0;
	ECC_st EC_P = { {0x00},0x00 };
	ECC_Jacobian_st EC_J_Temp= { {0x00},0x00 };
	bigint_st BN_Temp = { {0x00},0x00 };
	word P[ECC_256] = { 0x00, };

	int first_bit = 0;
	Set_integer(&BN_Temp, 0, P);
	Set_ECC(&EC_P, 0, P, P);
	Jacobian_To_Affine(EC_Base, &EC_Temp, Prime);
	NAFw_Fuction(Scalar, NAFw, Prime);
	Compute_Table_Jacobian(&EC_Temp, Table, Prime); // Table Complete


	for (cnt_i = 256; cnt_i >= 0; cnt_i--) {

		temp = (NAFw[cnt_i]);

		if (temp != 0) {
			first_bit += 1;
			if (first_bit == 1) {
				if (temp > 0) {
					ECC_copy(&EC_P, &Table[(temp - 1) / 2]);
					Affine_To_Jacobian(&EC_P, &EC_J_Temp, Prime);
					
				}
				else {
					ECC_copy(&EC_P, &Table[(-temp - 1) / 2]);
					Affine_To_Jacobian(&EC_P, &EC_J_Temp, Prime);
					
				}

			}

			else {
				if (temp > 0) {
				
					Jacobian_DBL(&EC_J_Temp, EC_R,Prime);
					Jacobian_ADD(EC_R, &Table[(temp - 1) / 2], &EC_J_Temp, Prime);
		
					
				}
				else {
			
					Jacobian_DBL(&EC_J_Temp, EC_R, Prime);
					ECC_copy(&EC_Temp, &Table[(-temp - 1) / 2]);
					Subtraction(Prime, &(EC_Temp.y), &BN_Temp);
					Subtraction_in_Fp(&BN_Temp, Prime, &(EC_Temp.y));
					Jacobian_ADD(EC_R, &EC_Temp, &EC_J_Temp, Prime);
			

			
				}

			}
		}
		else {
			if ((first_bit > 0)) {
				Jacobian_DBL(&EC_J_Temp, EC_R, Prime);

				Jacobian_Copy(&EC_J_Temp, EC_R);
			
			}

		}
		Jacobian_Copy(EC_R,&EC_J_Temp);
	}


	return;

}
void Comb_Table(ECC_st* EC_Base, ECC_st EC_Table[], unsigned char Table[][32], bigint_st* Scalar, bigint_st* Prime) { // A는 32칸을 가짐

	int temp = 0;
	int cnt_i = 0;
	int cnt_j = 0;
	int cnt_k = 0;
	int first_bit = 0;

	ECC_Jacobian_st EC_J_Temp = { {0x00},0x00 };
	ECC_Jacobian_st EC_C_Temp = { {0x00},0x00 };
	ECC_Jacobian_st Temp2 = { {0x00},0x00 };


	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		for (cnt_j = 31, cnt_k = 0; (cnt_j >= 0); cnt_j--, cnt_k++) {
			Table[cnt_i][cnt_k] = ((Scalar->a[cnt_i]) >> cnt_j) & 1;
		}
	} // 테이블이 만들어짐
	//Table 계산을 위한 값 필요함


	ECC_copy(&EC_Table[0], EC_Base);
	Affine_To_Jacobian(EC_Base, &EC_J_Temp, Prime);



	for (cnt_k = 1; cnt_k < 8; cnt_k++) {
		for (cnt_i = 0; cnt_i < 32; cnt_i++) {
			Jacobian_DBL(&EC_J_Temp, &EC_C_Temp, Prime);
			Jacobian_Copy(&EC_J_Temp, &EC_C_Temp);
		}


		Jacobian_To_Affine(&EC_C_Temp, &EC_Table[cnt_k], Prime);

	}// 테이블 끝.



	return;
}

void Comb_Scalar_Multiplication(ECC_Jacobian_st* EC_Base, ECC_Jacobian_st* EC_R,ECC_st EC_Table[], unsigned char Table[][32], bigint_st* Scalar, bigint_st* Prime) {
	
	int cnt_i = 0;
	int first_bit = 0;
	int cnt_j = 0;
	ECC_Jacobian_st EC_C_Temp= { {0x00},0x00 };
	ECC_Jacobian_st EC_J_Temp = { {0x00},0x00 };


	for (cnt_i = 0; cnt_i < 32; cnt_i++) {
	
		for (cnt_j = 0; cnt_j < ECC_256; cnt_j++) {
			if (Table[cnt_j][cnt_i] == 1) {
				first_bit += 1;
			//	printf("cnt_j=%d\n", cnt_j);
				if (first_bit == 1) {
					Affine_To_Jacobian(&EC_Table[cnt_j], &EC_J_Temp, Prime);
					
					Jacobian_Copy(&EC_C_Temp, &EC_J_Temp);
				}
				else {
				//	printf("첫번째 비트가 아닐때 1 cnt_j=%d\n", cnt_j);
				
					Jacobian_ADD(&EC_C_Temp, &EC_Table[cnt_j], &EC_J_Temp, Prime);
					Jacobian_Copy(&EC_C_Temp, &EC_J_Temp);

				//	printf("After\n");
				//	Show_ECC_Jacobian(&EC_C_Temp);
				
					
				}
			}

		}
	//	printf("\n\n");
		
		Jacobian_DBL(&EC_C_Temp, &EC_J_Temp, Prime);
		Jacobian_Copy(EC_R, &EC_C_Temp);
		Jacobian_Copy(&EC_C_Temp, &EC_J_Temp);

	}
	return;

}










