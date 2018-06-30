#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <algorithm>    // pow
#include <cmath> //isan
#include <iomanip>
#include <iostream>
#include <fstream>
#include <time.h>
#include <ctime>

using namespace std;

#define EPSILON 0.01
#define MAX_ITERS 40
#define SHOW_ITERS 0
#define CREATE_TASKFILE 0
#define NUMBER_OF_VARS_AND_MOD 10,30
#define WRITE_TO_FILE 0
#define VECTOR_SIZE_BUFFER 100
int details = 1;

#define DEFAULT_TASKFILE_PATH "task.txt"
#define CREATE_TASKFILE_PATH "created_task.txt"
#define OUTPUTFILE_PATH "test.txt"

void vector_out(vector<int> vec);
void vector_out(vector<bool> vec);
void vector_out(vector< vector<int> > vec);
void vector_out(vector< vector<bool> > vecvec);
void testfile(unsigned int size, unsigned int modulus);

unsigned tmp;


void create_LEsystems(
		vector<vector<int> > *constrains,
		vector<vector<bool> > *binary_pos,
		vector<vector<int> > *Gleichungssysteme);
void create_binary_possibilities(
		vector<int> *costfunction,
		int A,
		vector<vector<bool> > *binary_pos);
int calc_number_terminated_variables(vector<vector<int> > constrains, int variables_counter);
int count_zeros(vector<bool> bool_vector);
vector<bool> convert_int2boolvector(int x);
void jacobiCalcDisplay( int numUnKnowns, double** mat, double* result);
void showXcheck( int num, double** mat, double* v1, double* v2 );
double dotProd( int num, double* v1, double* v2 );
int checkFlags( int num, int flags[] );
int checkIfDD( int numUnKnowns, double** mat );
void gaus(double** mata, int variables_counter, double* result);

struct final_result{
	double result;
	vector<bool> binaries;
};

int main(int argc, char** argv) {

	// init
	srand(time(NULL));
	cout.setf(ios::fixed);
	
	clock_t start = clock();

	streambuf *psbuf;
	ofstream filestr;
	streambuf* old_buffer = std::cout.rdbuf ();

	cout<<"Prozess started."<<endl;
	if(WRITE_TO_FILE){
		cout<<"Writing into file: "<<OUTPUTFILE_PATH<<endl;
		filestr.open (OUTPUTFILE_PATH);
		psbuf = filestr.rdbuf();
		cout.rdbuf(psbuf);
	}
	// Main-Variables
	string taskfile_path = DEFAULT_TASKFILE_PATH;
    string line,cost,tmp_string,cost_gloabl;

    int tmp_int = 0;
    int vorzeichen = 1;
    int variables_counter = 0;
    int number_terminated_variables = 0;
	int n = 0;
	int safemode = 1;

	int sol_counter = 0;

	double temp_res = 0;

    double **mat = NULL;
    double *costfunction_new = NULL;
    double* result = NULL;

    vector<int> costfunction;
    vector<int> tmp_function;
    vector<vector<int> > constrains;
    vector<vector<bool> > binary_pos;
    vector<vector<int> > gleichungssysteme;

	final_result final_result;

	// Handling incomming parameters
	if(argc >= 2){
		taskfile_path = argv[1];
		if(argc >= 3){
			safemode = atoi(argv[2]);
			details = atoi(argv[3]);
		}
	}

	if(CREATE_TASKFILE){
		testfile(NUMBER_OF_VARS_AND_MOD);
		taskfile_path = CREATE_TASKFILE_PATH;
		cout<< float( clock() - start)/ CLOCKS_PER_SEC<<endl;
	}

	ifstream readFile (taskfile_path);
	cout<<"Using taskfile: "<<taskfile_path<<"\n"<<endl;

	if(details == 1)cout<<"Details are shown. Maybe cause a lot of output lines.\n"<<endl;
	if(safemode == 0)cout<<"Attention! Running in UNSAFE Mode! Diagonal dominance will not be checked."<<endl;


// Einlesen der Datei, Zeile für Zeile, Auswertung Zeichen für Zeichen
    while(getline(readFile,line))
    {
    	if(line[0] == '/' || line[0] == ' ') continue;

        if(line[0] == 'm'){
        	if(line[1] == 'i'){
        		cost = "min";
        	}
        	else if(line[1] == 'a'){
        		cost = "max";
        	}
        	cost_gloabl=cost;
        }

        for(unsigned i = 0; i <= line.size(); i++){
        	if(isdigit(line[i]) && !isalpha(line[i-1])) tmp_string = tmp_string + line[i];

        	else if( !tmp_string.empty() ) {
        		tmp_int = atoi(tmp_string.c_str());
        		tmp_int *= vorzeichen;
        		tmp_function.push_back(tmp_int);
        		tmp_int = 0;
        		tmp_string.clear();
        		vorzeichen = 1;
        		if(cost=="max" || cost=="min")variables_counter++;
        	}

        	else if( line[i] == '-') vorzeichen = -1;
        	else if( line[i] == '+') vorzeichen = 1;
        	else if( line[i] == '<' || line[i] =='>'){
        		for(std::vector<vector<int> >::const_iterator it = constrains.begin(); it != constrains.end(); it++){
        			tmp_function.push_back(0);
        		}
        		tmp_function.push_back(1);
        		costfunction.push_back(0);
        	}
        }

        if(cost == "max") {
        	costfunction = tmp_function;
        	cost.clear();
        }
        else if(cost == "min"){
        	costfunction = tmp_function;
        	cost.clear();
        }
        else if(!(tmp_function.empty())){
        	constrains.push_back(tmp_function);
        }
        tmp_function.clear();
    }

	cout<< float( clock() - start)/ CLOCKS_PER_SEC<<endl;


    if(cost_gloabl == "max")costfunction.push_back(1);
    else if (cost_gloabl == "min")costfunction.push_back(-1);

        vector< vector<int> >::iterator row;
        for (row = constrains.begin(); row != constrains.end(); row++) {
        	while((*row).size()<costfunction.size())(*row).insert(--(*row).end(),0);
        }

    number_terminated_variables = calc_number_terminated_variables(constrains,variables_counter);

    cout<<"Costfunction: "<<endl;
    vector_out(costfunction);
    cout<<"Constrains: "<<endl;
    vector_out(constrains);


    //--------------------------------------------------------------------------

    costfunction_new = (double*) malloc ((variables_counter+1)*sizeof(double));
    mat = (double**) malloc (variables_counter*sizeof(double*) );

    for(int u = 0; u < variables_counter; u++){
	mat[u] = (double*) malloc ( (variables_counter+1)*sizeof(double) );
    }

    result = (double*) malloc ( (variables_counter)*sizeof(double));

    // Berechnungschleife
    cout<<"\nStarting calculation iterations..."<<endl;
    while(tmp <= pow(2,costfunction.size()-1)){
    create_binary_possibilities(&costfunction, number_terminated_variables, &binary_pos);

    create_LEsystems(&constrains, &binary_pos, &gleichungssysteme);

    if(details==1 && binary_pos.size() !=0){
    	cout<<"Binary possibilities:"<<endl;
    	vector_out(binary_pos);
    	cout<<"Resulting equationsystems:"<<endl;
    	vector_out(gleichungssysteme);
    }

    for(vector<vector<bool> >::iterator bool_row_ite = binary_pos.begin(); bool_row_ite != binary_pos.end(); ++bool_row_ite){
    	for(unsigned bool_col = 0; bool_col <(*bool_row_ite).size(); ++bool_col)
    	{
    		if((*bool_row_ite).at(bool_col) != 0){
    			for(int i = 0; i < variables_counter; i++){
    				mat[i][n] = constrains[i].at(bool_col);
    				costfunction_new[n] = costfunction.at(bool_col);
    			}
				n++;
    		}
    	}
    	n=0;

    	if(safemode==1){
    		if(checkIfDD(variables_counter, mat)){
    			jacobiCalcDisplay(variables_counter, mat, result);
    		}
    		else gaus(mat, variables_counter, result);
		}
    	else {
    		jacobiCalcDisplay(variables_counter, mat, result);
    	}
		for(int i=0; i < variables_counter; i++){
				if(result[i]<0){
					temp_res = 0;
					if(details)cout<<"\nNot a valid solution."<<endl;
					break;
				}
				else if(result[i]>0){
				temp_res += costfunction_new[i]*result[i];
			}
		}

		if(temp_res > final_result.result){
			final_result.result = temp_res;
			final_result.binaries = *bool_row_ite;
		    cout<<"\nBest solution: "<<final_result.result<<endl;
		    vector_out(final_result.binaries);
		    sol_counter++;
		    cout<<endl;
		}
		else if(details) cout<<"\nSolution: "<<temp_res<<endl;
    }
    binary_pos.clear();
    gleichungssysteme.clear();
    }

    free(mat);

    if(sol_counter == 0)cout<<"No solution found."<<endl;
    else cout<<"Number of solutuons: "<<sol_counter<<endl;

    if(WRITE_TO_FILE){
    	filestr.close();
    	std::cout.rdbuf (old_buffer);
    }
	cout<<"Process took "<< float( clock() - start)/ CLOCKS_PER_SEC<<" seconds."<<endl;
    cout<<"\nCalculation finished."<<endl;

    return 0;
}

int checkIfDD( int numUnKnowns, double** mat )
{
    int   m, n, dd = 0;
    double* chkdd;

    double* sumdd = (double*) malloc( numUnKnowns*sizeof(double) );
    chkdd = (double*) malloc( numUnKnowns*sizeof(double) );
    for( m = 0 ; m < numUnKnowns ; m++ )
        chkdd[m] = sumdd[m] = 0; /* all set to zero ... */
    if(details)cout<<"Checking if the matrix is (strictly) diagonally dominant..."<<endl;
    for( m = 0 ; m < numUnKnowns ; m++ )
    {
        for( n = 0 ; n < numUnKnowns ; n++ )
        {
            sumdd[m] += fabs(mat[m][n] );
        }
        sumdd[m] -= fabs(mat[m][m]);
        chkdd[m] = fabs(mat[m][m]);
        if(chkdd[m] > sumdd[m])
        {
            if(details)cout<<chkdd[m]<<" >  "<<sumdd[m]<<endl;
            dd++;
        }
        else
        	if(details)cout<<chkdd[m]<<" !> "<<sumdd[m]<<endl;
    }

    if(dd == numUnKnowns)
    {
        if(details)cout<<"\nYes, the matrix is (strictly) diagonally dominant."<<endl;
    }
    else
    {
        if(details)cout<<"\nNo, the matrix is NOT (strictly) diagonally dominant ... using Gaus."<<endl;
        free( sumdd );
        free( chkdd );
        return 0; /* false */
    }
    free( sumdd );
    free( chkdd );
    return 1; /* true */
}


void jacobiCalcDisplay( int numUnKnowns, double** mat, double* res)
{
	bool nan = false;
	bool max_it_exceeded = false;
    int* flag;
    //double* res;
    //res = (double*) malloc( numUnKnowns*sizeof(double) );
    int i, j, counter = 0;
    double* var = (double*) malloc( numUnKnowns*sizeof(double) );
    flag = (int*) malloc( numUnKnowns*sizeof(int) );
    for(i = 0 ; i < numUnKnowns ; i++ )
        var[i] = res[i] = flag[i] = 0;

    if(details){
    	cout<<"\nSTART JACOBI CALCULATING"<<endl;
    	cout<<"The initial value of each array element was set to zero ..."<<endl;
    }
    do
    {
        counter++;
        /* for each iteration keep a copy of the old results ... */
        for(i = 0 ; i < numUnKnowns ; i++ )
        {
            var[i] = res[i];
        }
        if( SHOW_ITERS ) cout<<"\nIteration number "<<counter<<endl;
        for(i = 0 ; i < numUnKnowns ; i++ ) /* calculation */
		{
            res[i] = mat[i][numUnKnowns];
            for(j = 0 ; j < numUnKnowns ; j++ )
                res[i] = res[i] - mat[i][j]*var[j] ;
            res[i] = res[i] + mat[i][i]*var[i] ;
            res[i] = res[i] / mat[i][i] ;
            if( SHOW_ITERS ) cout<<"x"<<i<<" = "<<res[i]<<endl;
            if( fabs(res[i] - var[i]) < EPSILON) /* stop condition */
                flag[i]++;
            if( isnan(res[i])){
            	nan = true;
            	flag[i]++;

            }
            if( counter==MAX_ITERS) {
            	max_it_exceeded = true;
            	flag[i]++;
            }
        }
    }while( !checkFlags( numUnKnowns, flag ) );
    if(details)cout<<"\nThe RESULTS of "<<counter<<" ITERATIONS"<<endl;
    /*  cross check ...*/
    for( i = 0 ; i < numUnKnowns ; i++)
    {
        var[i] = dotProd( numUnKnowns, mat[i], res );
    }
    if(details)showXcheck( numUnKnowns, mat, res, var );

    if(nan){
    	//cout<<"\nMin one result is not-a-number..."<<endl;
    	nan = false;
    }

    if(max_it_exceeded){
    	//cout<<"\nMax iterations exceeded..."<<endl;
    	max_it_exceeded = false;
    }
    /* show sol'n vector (again) ... and free up all dynamic memory  */
    //cout<<"\nSolution vector"<<endl;
    for( i = 0 ; i < numUnKnowns ; i++)
    {
        //cout<<"x"<<i<<" = "<<res[i]<<endl;
        //free(mat[i]);
    }
    //result = res;
    //free( mat );
    free( flag );
    //free( res );
    free( var );
}


void showXcheck( int num, double** mat, double* v1, double* v2 )
{
    int i, j;
    cout<<"\nCross checking ... \nMatrix times sol'n vector ="<<
         " cal. vector vs. original RHS vector"<<endl;
    for( i = 0 ; i < num ; i++)
    {
        cout<<"|";
        for( j =0 ; j < num ; j++ )
            cout<<mat[i][j]<<" ";
        cout<<"| |"<<v1[i]<<"| |"<<v2[i]<<"| vs |"<<mat[i][num]<<endl;
    }
}

double dotProd( int num, double* v1, double* v2 )
{
    int i;
    double sum =0;
    for( i=0; i<num; ++i ) sum += v1[i]*v2[i];
    return sum;
}

int checkFlags( int num, int flags[] )
{
    int i;
    for( i=0; i<num; ++ i)
        if( flags[i] == 0 ) return 0;
    return 1;
}


void create_binary_possibilities(vector<int> *costfunction, int A, vector<vector<bool> > *binary_pos){

    vector<bool>  bool_vector;

    for( ; tmp <= pow(2,costfunction->size()-1); tmp++)
    {
    	bool_vector = convert_int2boolvector(tmp);

    	while(bool_vector.size()< costfunction->size()-1)bool_vector.push_back(0);

    	if(count_zeros(bool_vector)== A && binary_pos->size() <= VECTOR_SIZE_BUFFER){
    		bool_vector.push_back(1);
    		binary_pos->push_back(bool_vector);
    		if(binary_pos->size() == VECTOR_SIZE_BUFFER){
    			tmp++;
    			return;
    		}
    	}
    }
}

void create_LEsystems(
		vector<vector<int> > *constrains,
		vector<vector<bool> > *binary_pos,
		vector<vector<int> > *Gleichungssysteme){

	vector<int> tmp_vector;

    for(unsigned t = 0; t < binary_pos->size();t++){
    	for(vector<vector<int> >::iterator it = constrains->begin(); it != constrains->end(); ++it ){
    		tmp_vector = *it;
    		for(unsigned g = 0; g < constrains->at(0).size(); g++){
    			tmp_vector.at(g) = (*it).at(g) * binary_pos->at(t).at(g);
    		}
    		Gleichungssysteme->push_back(tmp_vector);
    		tmp_vector.clear();
    	}
    }
}

vector<bool> convert_int2boolvector(int x) {
  vector<bool> tmp_bool;
  while(x) {
    if (x&1)
      tmp_bool.push_back(1);
    else
      tmp_bool.push_back(0);
    x>>=1;
  }
  return tmp_bool;
}

int count_zeros(vector<bool> bool_vector)
{
	int tmp = 0;
	for(vector<bool>::iterator iti = bool_vector.begin(); iti != bool_vector.end(); iti++)
	{
		if((*iti) == 0) tmp++;
	}
	return tmp;
}


int calc_number_terminated_variables(vector<vector<int> > constrains, int variables_counter)
{
	int AH = constrains[0].size()-(variables_counter+1); //number of slack-variables
	int AG = constrains.size(); // number of equations
	int A = variables_counter+AH-AG; // number of null-variables
	cout<<"Variables: " <<variables_counter<<" Slacks: " << AH <<" Equis: "<<AG<<" 0-Variables: "<<A<<endl;
	return A;
}

void vector_out(vector<int> vec)
{
	for (std::vector<int>::const_iterator i = vec.begin(); i != vec.end(); ++i)
	    std::cout << *i << ' ';
	cout<<endl;
}

void vector_out(vector<bool> vec)
{
	for (vector<bool>::const_iterator i = vec.begin(); i != vec.end(); ++i)cout << *i << ' ';
	cout<<endl;
}

void vector_out(vector< vector<int> > vecvec)
{
	for (vector< vector<int> >::iterator row = vecvec.begin(); row != vecvec.end(); row++)vector_out(*row);
	cout<<endl;
}

void vector_out(vector< vector<bool> > vecvec)
{
	for (vector< vector<bool> >::iterator row = vecvec.begin(); row != vecvec.end(); row++) vector_out(*row);
	cout<<endl;
}

void gaus(double** mata, int variables_counter, double* x){
	int n = variables_counter;
	int i,j,k = 0;

	double mat[n][n+1];

	for(int l = 0; l < variables_counter; l++){
		for(int m = 0; m<=variables_counter;m++){
			mat[l][m] = mata[l][m];
		}
	}

    if(details)cout<<"\nSTART GAUS ELIMINATION"<<endl;


	    for (i=0;i<n;i++)                    //Pivotisation
	        for (k=i+1;k<n;k++)
	            if (mat[i][i]<mat[k][i])
	                for (j=0;j<=n;j++)
	                {
	                    double temp=mat[i][j];
	                    mat[i][j]=mat[k][j];
	                    mat[k][j]=temp;
	                }

	    if(details){
			cout<<"\nThe matrix after Pivotisation is:\n";
			for (i=0;i<n;i++)            //print the new matrix
			{
				for (j=0;j<=n;j++){

					cout<<mat[i][j]<<" ";
					if(j==n-1)cout<<" | ";
				}
				cout<<"\n";
			}
	    }

	    for (i=0;i<n-1;i++)            //loop to perform the gauss elimination
	        for (k=i+1;k<n;k++)
	            {
	                double t=mat[k][i]/mat[i][i];
	                for (j=0;j<=n;j++)
	                    mat[k][j]=mat[k][j]-t*mat[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
	            }

	    if(details){
			cout<<"\nThe matrix after gauss-elimination is as follows:\n";
			for (i=0;i<n;i++)            //print the new matrix
			{
				for (j=0;j<=n;j++){
					cout<<mat[i][j]<<" ";
					if(j==n-1)cout<<" | ";
				}
				cout<<"\n";
			}
	    }
	    for (i=n-1;i>=0;i--)                //back-substitution
	    {                        //x is an array whose values correspond to the values of x,y,z..
	        x[i]=mat[i][n];                //make the variable to be calculated equal to the rhs of the last equation
	        for (j=0;j<n;j++)
	            if (j!=i)            //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
	                x[i]=x[i]-mat[i][j]*x[j];
	        x[i]=x[i]/mat[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
	    }

	    if(details){
	    	cout<<"\nThe values of the variables are as follows:\n";
	    	for (i=0;i<n;i++)
	    		cout<<"X"<<i<<" = "<<x[i]<<endl; // Print the values of x, y,z,....
	    	cout<<endl;
	    }
}


void testfile(unsigned int size, unsigned int modulus)
{
    cout << "Start generating Taskfile." << endl;
    if(modulus == 0 || modulus == 1) modulus = 3;
    ofstream writing_stream;
    writing_stream.open(CREATE_TASKFILE_PATH);

    writing_stream << "max: " << (rand()%modulus) << "*x1";
    for(unsigned int a = 2; a <= size; a++)
    {
        if(rand()%2 == 0) writing_stream << " + " << (rand()%modulus) << "*x" << a%9 ;
        else writing_stream << " + " << (rand()%modulus) << "*x" << a%9 ;
    }
    writing_stream <<";"<< endl;

    for(unsigned int i = 1; i <= size; i++)
    {
        for(unsigned int j = 1; j <= size; j++)
        {
            if(i == j)
            {
                if(j == 1)
                    {
                        if(rand()%2 == 0) writing_stream << (rand()%modulus) << "*x1";
                        else writing_stream << "" << (rand()%modulus) << "*x1";
                    }
                else
                    {
                        if(rand()%2 == 0) writing_stream << " + " << (rand()%modulus) << "*x" << j%9 ;
                        else writing_stream << " - " << (rand()%modulus) << "*x" << j%9 ;
                    }
            }
            else
            {
                if(j == 1)
                    {
                        if(rand()%2 == 0) writing_stream << rand()%(modulus) << "*x1";
                        else writing_stream << "" << rand()%(modulus) << "*x1";
                    }
                else
                {
                    if(rand()%2 == 0) writing_stream << " + " << rand()%(modulus) << "*x" << j%9 ;
                    else writing_stream << " - " << rand()%(modulus) << "*x" << j%9 ;
                }
            }
        }
        writing_stream << " <= " << rand()%modulus << ";\n";
    }
    cout << "Taskfile generated." << endl;
}
