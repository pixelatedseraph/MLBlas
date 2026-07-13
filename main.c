#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdint.h>
#include<stdbool.h>
#include<unistd.h>
#include<time.h>
#include<limits.h>

typedef enum{
    Double,
    Int,
}MatrixType;

typedef struct{
    void*       data;
    size_t      rows;
    size_t      cols;
    MatrixType  type;
}Matrix;

typedef enum {
    HeapAllocationError,
    InvalidDimensionError,
    TypeMismatchError,
}MatrixErrors;

typedef struct{
    bool            ok;
    MatrixErrors    matrix_error;
    Matrix          matrix;
}MatrixResult;


static inline void thread_panic(const char* fmt){
    fprintf(stderr,"[Linxel Error]: %s\n",fmt);
    exit(EXIT_FAILURE);
}

#define for_each(idx,matrix) for(size_t idx = 0; idx < matrix.rows * matrix.cols ; ++idx)

static inline MatrixResult matrix_create(MatrixType matrix_type,size_t rows,size_t cols){
    Matrix matrix;
    MatrixResult result;
    size_t matrix_size = 0;
    if(matrix_type == Double){
        matrix_size = sizeof(double);
    }
    else if(matrix_type == Int){
        matrix_size = sizeof(int);
    }
    else{
        matrix_size = 0;
        thread_panic("UnsupportedTypeError");
    }
    if(!(matrix.data = malloc(matrix_size * rows * cols))){
        result.ok = false;
        result.matrix_error = HeapAllocationError;
        return result;
    }   
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.type = matrix_type;

    result.ok = true;
    result.matrix = matrix;

    return result;
}

static inline void matrix_fill_int(Matrix matrix,int value){
    int* restrict matrix_data_int = (int*)matrix.data;

    for_each(i,matrix){
        matrix_data_int[i] = value;
    }

}

static inline void matrix_fill_random_int(Matrix matrix,int low,int high){
    int* restrict matrix_data_int = (int*)matrix.data;

    uint64_t range = (uint64_t)(int64_t)high - (uint64_t)(int64_t)low + 1;

    for_each(i,matrix){
        matrix_data_int[i] =  low + rand() % (range);
    }
}


static inline void matrix_fill_double(Matrix matrix,double value){
    double* restrict matrix_data_double = (double*)matrix.data;

    for_each(i,matrix){
        matrix_data_double[i] = value;
    }

}

static inline void matrix_fill_random_double(Matrix matrix,double low,double high){
    double* matrix_data_double = (double*)matrix.data;

    for_each(i,matrix){
        double double_value = (double)rand() / RAND_MAX;
        double filler_value = low + (high - low) * double_value;

        matrix_data_double[i] = filler_value;
    }
}

static inline size_t get_linear_index(Matrix matrix,size_t row,size_t col){
    if (matrix.rows <= row || matrix.cols <= col) {
        thread_panic("InvalidDimensionsError");
    }
    return row * matrix.cols + col;
}


static inline void matrix_print(const Matrix matrix){
    if(matrix.type == Int){
        int* restrict matrix_data_int = (int*)matrix.data;

        for(size_t i = 0 ; i < matrix.rows ; ++i){
            for(size_t j = 0 ; j < matrix.cols ; ++j){
                printf("%d ",matrix_data_int[get_linear_index(matrix,i,j)]);
            }
            printf("\n");
        }
        printf("\n");
    }
    else if(matrix.type == Double){
        double* restrict matrix_data_double = (double*)matrix.data;

        for(size_t i = 0 ; i < matrix.rows ; ++i){
            for(size_t j = 0 ; j < matrix.cols ; ++j){
                printf("%.2f ",matrix_data_double[get_linear_index(matrix,i,j)]);
            }
            printf("\n");
        }
        printf("\n");
    }
    else{
        thread_panic("UnprintableMatrixError");
    }

}

static inline MatrixResult matrix_add(const Matrix lhs,const Matrix rhs){
    MatrixResult result;
    Matrix       matrix;

    if((lhs.rows != rhs.rows) || (lhs.cols != rhs.cols)){
        result.ok = false;
        result.matrix_error = InvalidDimensionError;
        return result;
    }

    if(lhs.type != rhs.type){
        result.ok = false;
        result.matrix_error = TypeMismatchError;
        return result;
    }

    if(lhs.type == Int){
        int* restrict lhs_matrix_data_int = (int*)lhs.data;
        int* restrict rhs_matrix_data_int = (int*)rhs.data;

        MatrixResult internal_result = matrix_create(Int,lhs.rows,lhs.cols);

        if(internal_result.ok == false){
            return internal_result;
        }

        matrix = internal_result.matrix;
        int* restrict result_matrix_data_int = (int*)matrix.data;

        for(size_t i = 0 ; i < lhs.rows * lhs.cols ; ++i){
            result_matrix_data_int[i] = lhs_matrix_data_int[i] + rhs_matrix_data_int[i];            
        }

        result = internal_result;
        return result;
    }
    else{
        double* restrict lhs_matrix_data_double = (double*)lhs.data;
        double* restrict rhs_matrix_data_double = (double*)rhs.data;

        MatrixResult internal_result = matrix_create(Double,lhs.rows,lhs.cols);

        if(internal_result.ok == false){
            return internal_result;
        }

        matrix = internal_result.matrix;
        double* restrict result_matrix_data_double = (double*)matrix.data;

        for(size_t i = 0 ; i < lhs.rows * lhs.cols ; ++i){
            result_matrix_data_double[i] = lhs_matrix_data_double[i] + rhs_matrix_data_double[i];            
        }

        result = internal_result;
        return result;
    }
}

static inline MatrixResult matrix_sub(const Matrix lhs,const Matrix rhs){
    MatrixResult result;
    Matrix       matrix;

    if((lhs.rows != rhs.rows) || (lhs.cols != rhs.cols)){
        result.ok = false;
        result.matrix_error = InvalidDimensionError;
        return result;
    }

    if(lhs.type != rhs.type){
        result.ok = false;
        result.matrix_error = TypeMismatchError;
        return result;
    }

    if(lhs.type == Int){
        int* restrict lhs_matrix_data_int = (int*)lhs.data;
        int* restrict rhs_matrix_data_int = (int*)rhs.data;

        MatrixResult internal_result = matrix_create(Int,lhs.rows,lhs.cols);

        if(internal_result.ok == false){
            return internal_result;
        }

        matrix = internal_result.matrix;
        int* restrict result_matrix_data_int = (int*)matrix.data;

        for(size_t i = 0 ; i < lhs.rows * lhs.cols ; ++i){
            result_matrix_data_int[i] = lhs_matrix_data_int[i] - rhs_matrix_data_int[i];            
        }

        result = internal_result;
        return result;
    }
    else{
        double* restrict lhs_matrix_data_double = (double*)lhs.data;
        double* restrict rhs_matrix_data_double = (double*)rhs.data;

        MatrixResult internal_result = matrix_create(Double,lhs.rows,lhs.cols);

        if(internal_result.ok == false){
            return internal_result;
        }

        matrix = internal_result.matrix;
        double* restrict result_matrix_data_double = (double*)matrix.data;

        for(size_t i = 0 ; i < lhs.rows * lhs.cols ; ++i){
            result_matrix_data_double[i] = lhs_matrix_data_double[i] - rhs_matrix_data_double[i];            
        }

        result = internal_result;
        return result;
    }
}

static inline MatrixResult matrix_clone(const Matrix matrix){
    MatrixResult result;

    result = matrix_create(matrix.type,matrix.rows,matrix.cols); 
    if(result.ok == false){
        return result;
    }

    result.matrix.rows = matrix.rows;
    result.matrix.cols = matrix.cols;
    result.matrix.type = matrix.type;

    if(matrix.type == Int){
        int* restrict matrix_data_int = (int*)matrix.data;
        int* restrict result_matrix_int = (int*)result.matrix.data;
        
        for_each(i,matrix){
            result_matrix_int[i] = matrix_data_int[i];
        }
        return result;
    }
    else{
        double* restrict matrix_data_double = (double*)matrix.data;
        double* restrict result_matrix_double = (double*)result.matrix.data;
        
        for_each(i,matrix){
            result_matrix_double[i] = matrix_data_double[i];
        }
        return result; 
    }
}

static inline bool matrix_equal(const Matrix lhs, const Matrix rhs){
    if(lhs.type != rhs.type)
        return false;
    if((lhs.rows != rhs.rows) || (lhs.cols != rhs.cols))
        return false;
    
    if(lhs.type == Int){
        int* restrict lhs_data_int = (int*)lhs.data;
        int* restrict rhs_data_int = (int*)rhs.data;
        
        for(size_t i = 0 ; i < lhs.rows * lhs.cols ; ++i){
            if(lhs_data_int[i] != rhs_data_int[i]){
                return false;
            }
        }
        return true;
    }
    else{
        double* restrict lhs_data_int = (double*)lhs.data;
        double* restrict rhs_data_int = (double*)rhs.data;
        
        for(size_t i = 0 ; i < lhs.rows * lhs.cols ; ++i){
            if(lhs_data_int[i] != rhs_data_int[i]){
                return false;
            }
        }
        return true;
    }
}

static inline void matrix_cleanup(Matrix* restrict matrix){
    if(matrix){
        matrix->rows = 0;
        matrix->cols = 0;
        free(matrix->data);
        matrix->data = NULL;
    }
}



#define drop __attribute__((cleanup(matrix_cleanup))) // automatically calls free when goes out of scope

int main(int argc, char const *argv[]){
    (void)argc; (void) argv;
    srand(time(NULL));

   

    for(int i = 0 ; i < 100000000 ; ++i){
        drop Matrix m1 =  matrix_create(Int,10,4).matrix; // sort of like unwrap() <- use with extreme caution and never in production
        matrix_fill_random_int(m1,INT_MIN,INT_MAX);
        drop Matrix m2 =  matrix_clone(m1).matrix;

        if(matrix_equal(m1,m2) == false){
            matrix_print(m1);
            matrix_print(m2);
            fprintf(stderr,"Check for yourself buddy\n");
            return 69;
        }
    }

    return 0;
}
 


