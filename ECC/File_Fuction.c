#include"ECC.h"

void P_Read(FILE* fp, bigint_st* P) {
	int temp = 0x00;
	int cnt_i = 0;

	for (cnt_i = ECC_256-1; cnt_i >= 0; cnt_i--) {
		fscanf(fp, "%08X", &temp);
		P->a[cnt_i] = temp;
	}

	return;
}
void P_Read_Mul(FILE* fp, bigint_st_mul* P) {
	int temp = 0x00;
	int cnt_i = 0;

	for (cnt_i = ECC_256+ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		fscanf(fp, "%08X", &temp);
		P->a[cnt_i] = temp;
	}

	return;
}

void File_Read(FILE* fp, FILE* fp2, bigint_st* A, bigint_st* B) {
	int cnt_i = 0;
	word temp[ECC_256] = { 0x00, };
	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		fscanf(fp, "%08X", &temp[cnt_i]);
		A->a[cnt_i] = temp[cnt_i];
	}
	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		fscanf(fp2, "%08X", &temp[cnt_i]);
		B->a[cnt_i] = temp[cnt_i];
	}
	return;
}

void File_Write(FILE* fp, bigint_st* NUMBER) {
	int cnt_i = 0;
	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		fprintf(fp, "%08X", NUMBER->a[cnt_i]);
	}	fprintf(fp, "\n\n");
	return;

}


void File_Write_ECC(FILE* fp, ECC_st* NUMBER) {
	int cnt_i = 0;
	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		fprintf(fp, "%08X", NUMBER->x.a[cnt_i]);
	}	fprintf(fp, "\n");
	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		fprintf(fp, "%08X", NUMBER->y.a[cnt_i]);
	}	fprintf(fp, "\n\n");
	return;

}
void File_Write_Jacobian(FILE* fp, ECC_Jacobian_st* NUMBER) {
	int cnt_i = 0;
	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		fprintf(fp, "%08X", NUMBER->x.a[cnt_i]);
	}	fprintf(fp, "\n");
	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		fprintf(fp, "%08X", NUMBER->y.a[cnt_i]);
	}	fprintf(fp, "\n");
	for (cnt_i = ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		fprintf(fp, "%08X", NUMBER->z.a[cnt_i]);
	}	fprintf(fp, "\n\n");
	return;

}

void File_Write_Mul(FILE* fp, bigint_st_mul* NUMBER) {
	int cnt_i = 0;
	for (cnt_i = ECC_256+ECC_256 - 1; cnt_i >= 0; cnt_i--) {
		fprintf(fp, "%08X", NUMBER->a[cnt_i]);
	}	fprintf(fp, "\n\n");
	return;

}