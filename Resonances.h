
#define MAX_ORDER 3

typedef struct ResonanceData {
	int Nres,MaxOrder;
	int* ResonanceIndices;
	double* ResonanceCoefficients;
} ResonanceData;


void initialize_ResonanceData( ResonanceData * r);
void AddResonance(ResonanceData* r, int res_j, int order, int epower,double a);
void free_ResonanceData( ResonanceData * r);