#define true 1
#define false 0
#define MAX_ORDER 8
#define MAX_J 100
void initialize_ResonanceData( ResonanceData * r, int include0th);
void AddResonance(ResonanceData* r, int res_j, int order, int epower,double a);
void free_ResonanceData( ResonanceData * r);
