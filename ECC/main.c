#include "ECC.h"

int64_t cpucycles(void) { return __rdtsc(); }
#if Controler==0
int main() {
	int cnt_i = 0;
	int cnt_j = 0;
	int compare = 0;
	unsigned long long cycles, cycles1, cycles2;
	unsigned long long total = 0x00;


	FILE* P_rfp = NULL;
	FILE* A_rfp = NULL;
	FILE* B_rfp = NULL;
	FILE* rfp = NULL;
	FILE* Out_wfp = NULL;
	word* Mul_Result_p = NULL;
	word* P256 = NULL;
	bigint_st BigNumber = { {0x00},0x00 };
	bigint_st BigNumber2 = { {0x00},0x00 };
	bigint_st Mul_Result = { {0x00},0x00 };
	bigint_st Result = { {0x00},0x00 };
	bigint_st P = { {0x00},0x00 };
	bigint_st out = { {0x00},0x00 };

	Mul_Result_p = (word*)calloc(2 * ECC_256, sizeof(word));
	P256 = (word*)calloc(ECC_256, sizeof(word));


	fopen_s(&P_rfp, "P256��.txt", "rb");
	fopen_s(&A_rfp, "TV_opA_Inverse.txt", "rb");

	fopen_s(&rfp, "TV_PFINV_TV.txt", "rb");
	fopen_s(&Out_wfp, "TV_PFINV_TV(��갪Binary).txt", "wb");

	if (P_rfp == NULL) {
		printf("������ ������ �����ϴ�.");
	}

	//P�� �о��
	Set_integer(&P, 0, P256);

	P_Read(P_rfp, &P);




	while (cnt_i < 10000) {
		cnt_i++;
		//printf("%d��°\n", cnt_i);
		Set_integer(&BigNumber, 0,P256);

		Set_integer(&Result, 0, P256);
		Set_integer(&out, 0,  P256);

		P_Read(A_rfp, &BigNumber);

		//Show_integer(&BigNumber);

		cycles1 = cpucycles();
		Binary_Inversion(&BigNumber, &Result, &P);
		cycles2 = cpucycles();

		cycles = cycles2 - cycles1;
		total += cycles;
		P_Read(rfp, &out);
		//printf("��갪 : ");
		//Show_integer(&Result);
		//printf("������ : ");
		//Show_integer(&out);
	

		compare += Compare(&Result, &out);
		//printf("cycles: %10lld \n", cycles);

		File_Write(Out_wfp, &Result);

	}

	total = total / 10000;
	printf("��� ����Ŭ : %d\n", total);

	printf("���� �� ���� : %d\n", compare);
	fclose(P_rfp);
	fclose(A_rfp);

	fclose(rfp);
	fclose(Out_wfp);




	return 0;
}
#elif Controler==1
int main() {
	int cnt_i = 0;
	int cnt_j = 0;
	int compare = 0;
	unsigned long long cycles, cycles1, cycles2;
	unsigned long long total = 0x00;

	FILE* P_rfp = NULL;
	FILE* A_rfp = NULL;
	FILE* B_rfp = NULL;
	FILE* rfp = NULL;
	FILE* Out_wfp = NULL;
	FILE* Base_rfp = NULL;
	word* Mul_Result_p = NULL;
	word* P256 = NULL;
	bigint_st BigNumber = { {0x00},0x00 };
	bigint_st BigNumber2 = { {0x00},0x00 };
	bigint_st Mul_Result = { {0x00},0x00 };
	bigint_st Result = { {0x00},0x00 };
	bigint_st P = { {0x00},0x00 };
	bigint_st out = { {0x00},0x00 };
	bigint_st Scalar_Integer = { {0x00},0x00 };
	ECC_st EC_P = { {0x00},0x00 };
	ECC_st EC_Q = { {0x00},0x00 };
	ECC_st EC_R = { {0x00},0x00 };
	ECC_st EC_Base = { {0x00},0x00 };
	ECC_st temp = { {0x00},0x00 };

	Mul_Result_p = (word*)calloc(2 * ECC_256, sizeof(word));
	P256 = (word*)calloc(ECC_256, sizeof(word));

	fopen_s(&Base_rfp, "Base_Point.txt", "rb");
	fopen_s(&P_rfp, "P256��.txt", "rb");
	fopen_s(&A_rfp, "TV_SM.txt", "rb");

	fopen_s(&rfp, "TV_Scalar.txt", "rb"); // ����� 
	fopen_s(&Out_wfp, "TV_Scalar(��갪).txt", "wb"); // ��갪

	//P�� �о��
	Set_integer(&P, 0, P256);
	Set_ECC(&EC_Base, 0, P256, P256);
	P_Read(P_rfp, &P);
	//������ �б�
	P_Read(Base_rfp, &(EC_Base.x));
	P_Read(Base_rfp, &(EC_Base.y));
	Set_integer(&Scalar_Integer, 0, P256);

	
	while (cnt_i < 10000) {
		cnt_i++;
		printf("cnt_i=%d\n", cnt_i);
		Set_ECC(&EC_R, 0, P256, P256);
		Set_ECC(&EC_P, 0, P256, P256);
		Set_ECC(&EC_Q, 0, P256, P256);
		P_Read(rfp, &Scalar_Integer);

		//printf("������\n");
		//Show_ECC(&EC_Base);
		//printf("��Į��\n");
		//Show_integer(&Scalar_Integer);

		ECC_copy(&EC_P, &EC_Base);


		cycles1 = cpucycles();

		//Elliptic_Curve_DBL(&EC_P,&EC_Q, &EC_Base, &EC_R, &P);
		Left_to_Right(&EC_Base, &EC_R, &Scalar_Integer, &P);
		//Right_to_Left(&EC_Base, &EC_R, &Scalar_Integer, &P);
		cycles2 = cpucycles();


		cycles = cycles2 - cycles1;
		total += cycles;

		//	printf("�����:");
		//	Show_ECC(&EC_R);

		//	compare += Compare(Result, out);
			//printf("cycles: %10lld \n", cycles);

		File_Write_ECC(Out_wfp, &EC_R);
	}





	total = total / 100;
	printf("��� ����Ŭ : %d\n", total);
	printf("���� �� ���� : %d\n", compare);
	fclose(P_rfp);
	fclose(A_rfp);

	fclose(rfp);
	fclose(Out_wfp);




	return 0;
}

#elif Controler==2
int main() {
	int cnt_i = 0;
	int cnt_j = 0;
	int compare = 0;
	unsigned long long cycles, cycles1, cycles2;
	unsigned long long total = 0x00;

	FILE* P_rfp = NULL;
	FILE* A_rfp = NULL;
	FILE* B_rfp = NULL;
	FILE* rfp = NULL;
	FILE* Out_wfp = NULL;
	FILE* Base_rfp = NULL;
	word* Mul_Result_p = NULL;
	word* P256 = NULL;
	bigint_st BigNumber = { {0x00},0x00 };
	bigint_st BigNumber2 = { {0x00},0x00 };
	bigint_st Mul_Result = { {0x00},0x00 };
	bigint_st Result = { {0x00},0x00 };
	bigint_st P = { {0x00},0x00 };
	bigint_st out = { {0x00},0x00 };
	bigint_st Scalar_Integer = { {0x00},0x00 };
	ECC_st EC_P = { {0x00},0x00 };
	ECC_st EC_Q = { {0x00},0x00 };
	ECC_st EC_R = { {0x00},0x00 };
	ECC_st EC_Base = { {0x00},0x00 };
	ECC_st temp = { {0x00},0x00 };
	////////test code////////////
	ECC_Jacobian_st Test= { {0x00},0x00 };
	ECC_Jacobian_st Test2 = { {0x00},0x00 };

	////////////////////////////

	Mul_Result_p = (word*)calloc(2 * ECC_256, sizeof(word));
	P256 = (word*)calloc(ECC_256, sizeof(word));

	fopen_s(&Base_rfp, "Base_Point.txt", "rb");
	fopen_s(&P_rfp, "P256��.txt", "rb");
	fopen_s(&A_rfp, "TV_SM.txt", "rb");

	fopen_s(&rfp, "TV_Scalar.txt", "rb"); // ����� 
	fopen_s(&Out_wfp, "TV_Scalar(��갪).txt", "wb"); // ��갪

	//P�� �о��
	Set_integer(&P, 0, P256);
	Set_ECC(&EC_Base, 0, P256, P256);
	P_Read(P_rfp, &P);
	//������ �б�
	P_Read(Base_rfp, &(EC_Base.x));
	P_Read(Base_rfp, &(EC_Base.y));
	Set_integer(&Scalar_Integer, 0, P256);

	

	while (cnt_i < 10000) {
		cnt_i++;
		printf("Cnt_i=%d\n", cnt_i);
		Set_ECC(&EC_R, 0, P256, P256);
		Set_ECC(&EC_P, 0, P256, P256);
		Set_ECC(&EC_Q, 0, P256, P256);
		P_Read(rfp, &Scalar_Integer);

		//printf("������\n");
		//Show_ECC(&EC_Base);
		//printf("��Į��\n");
		//Show_integer(&Scalar_Integer);

	

		Affine_To_Jacobian(&EC_Base, &Test, &P);
	
		
		cycles1 = cpucycles();
		Jacobian_Left_to_Right(&Test, &Test2, &Scalar_Integer, &P);
		//Elliptic_Curve_Add(&EC_P,&EC_Q, &EC_Base, &EC_R, &P);
		//Left_to_Right(&EC_Base, &EC_R, &Scalar_Integer, &P);
		//Right_to_Left(&EC_Base, &EC_R, &Scalar_Integer, &P);
		//Jacobian_DBL(&Test, &Test2, &P);
		//Jacobian_ADD(&Test, &EC_Base, &Test2, &P);
		cycles2 = cpucycles();
		
		//Show_ECC_Jacobian(&Test2);
		Jacobian_To_Affine(&Test2, &EC_P, &P);

		cycles = cycles2 - cycles1;
		total += cycles;

			//printf("�����\n");
			//Show_ECC(&EC_P);

		//	compare += Compare(Result, out);
			//printf("cycles: %10lld \n", cycles);

		File_Write_ECC(Out_wfp, &EC_P);
	}





	total = total / 10000;
	printf("��� ����Ŭ : %lld\n", total);
	printf("���� �� ���� : %d\n", compare);
	fclose(P_rfp);
	fclose(A_rfp);

	fclose(rfp);
	fclose(Out_wfp);




	return 0;
}

#elif Controler==3
int main() {

	unsigned long long cycles, cycles1, cycles2;
	unsigned long long total = 0x00;


	bigint_st P = { {0x00},0x00 };
	bigint_st Scalar = { {0x00},0x00 };
	ECC_st EC_Base = { {0x00},0x00 };
	ECC_st EC_R = { {0x00},0x00 };
	FILE* fp = NULL;
	FILE* Base_rfp = NULL;
	FILE* Scalar_rfp = NULL;
	FILE* Out_wfp = NULL;
	word Zero[ECC_256] = { 0x00, };
	
	ECC_Jacobian_st Test = { {0x00},0x00 };
	ECC_Jacobian_st Test2 = { {0x00},0x00 };
	int cnt_i = 0;

	Set_ECC(&EC_Base, 0, Zero, Zero);
	Set_ECC(&EC_R, 0, Zero, Zero);
	Set_integer(&Scalar, 0, Zero);
	Set_integer(&P, 0, Zero);


	fopen_s(&fp, "P256��.txt", "rb");
	fopen_s(&Base_rfp, "Base_Point.txt", "rb");
	fopen_s(&Scalar_rfp, "TV_Scalar.txt", "rb");
	fopen_s(&Out_wfp, "TV_Scalar(��갪mNAF_Jacobian).txt", "wb"); // ��갪



	P_Read(fp, &P);

	P_Read(Base_rfp, &(EC_Base.x));
	P_Read(Base_rfp, &(EC_Base.y));

	while (cnt_i < 10000) {
		cnt_i++;
		//printf("%d\n", cnt_i);
		P_Read(Scalar_rfp, &Scalar);

		Affine_To_Jacobian(&EC_Base, &Test, &P);

		//Show_integer(&Scalar);
		
		cycles1 = cpucycles();
		NAF_Left_To_Right_Jacobian(&Test, &Test2, &Scalar, &P);
		//NAF_Left_To_Right_Affine(&EC_Base, &EC_R, &Scalar, &P);
		cycles2 = cpucycles();
		cycles = cycles2 - cycles1;
		total += cycles;
		Jacobian_To_Affine(&Test2, &EC_R, &P);
		//	Show_ECC(&EC_R);
	//	File_Write_ECC(Out_wfp, &EC_R);
	}

	total = total / 10000;
	printf("��� ����Ŭ : %lld\n", total);
}


#elif Controler==4
int main() {
	int cnt_i = 0;
	bigint_st P = { {0x00},0x00 };
	bigint_st Scalar = { {0x00},0x00 };
	ECC_st EC_Base = { {0x00},0x00 };
	ECC_st EC_R = { {0x00},0x00 };
	FILE* fp = NULL;
	FILE* Base_rfp = NULL;
	FILE* Scalar_rfp = NULL;
	FILE* Out_wfp = NULL;

	unsigned long long cycles, cycles1, cycles2;
	unsigned long long total = 0x00;

	word Zero[ECC_256] = { 0x00, };
	unsigned char Comb[ECC_256][32] = { 0x00, };

	ECC_Jacobian_st Test = { {0x00},0x00 };
	ECC_Jacobian_st Test2 = { {0x00},0x00 };
	ECC_st Table[8] = { {0x00},0x00 };

	Set_ECC(&EC_Base, 0, Zero, Zero);
	Set_ECC(&EC_R, 0, Zero, Zero);
	Set_integer(&Scalar, 0, Zero);

	Set_integer(&P, 0, Zero);


	fopen_s(&fp, "P256��.txt", "rb");
	fopen_s(&Base_rfp, "Base_Point.txt", "rb");
	fopen_s(&Scalar_rfp, "TV_Scalar.txt", "rb");
	//fopen_s(&Out_wfp, "TV_Scalar(��갪_Comb).txt", "wb"); // ��갪


	P_Read(fp, &P);
	P_Read(Base_rfp, &(EC_Base.x));
	P_Read(Base_rfp, &(EC_Base.y));

	while (cnt_i < 10000) {
		cnt_i++;
		printf("%d\n", cnt_i);
		P_Read(Scalar_rfp, &Scalar);


		Comb_Table(&EC_Base, Table,Comb, &Scalar, &P);

		Affine_To_Jacobian(&EC_Base, &Test, &P);
		cycles1 = cpucycles();
		Comb_Scalar_Multiplication(&Test, &Test2, Table,Comb, &Scalar, &P);
		
		cycles2 = cpucycles();

		cycles = cycles2 - cycles1;
		total += cycles;

		Jacobian_To_Affine(&Test2, &EC_R, &P);
		Show_ECC(&EC_R);
	//	File_Write_ECC(Out_wfp, &EC_R);
	}

	total = total / 10000;
	printf("��� ����Ŭ : %lld\n", total);
	
	fclose(fp);
	fclose(Base_rfp);
//	fclose(Out_wfp);
	fclose(Scalar_rfp);
	
	return 0;
}

#endif



/// ������� >>1 ĭ ¦Ȧ�� �������� 1���� Ȯ��
//��ä�κм��� �����ϰ��ϱ� ���� x�� ���� ���Ҷ� ������ r�� �����ְ� rx�� ������ ���� �� r�� �ٽ� ���ϸ� x�� ������ ���� ���ǿ�
//���� �ٸ��ϱ� r�� ���� �� �� ����.
//������ǥ�� ���� ��������� ����
// ������ ������ ������ �������� ��ü�ϰ����. -> �翵��ǥ��
//��ȣ �ٲܶ��� P-Y �������. 