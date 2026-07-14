#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdint.h>
#include<stdbool.h>
#include<unistd.h>
#include<time.h>
#include<limits.h>

//TODO 1 implement elipson equality 
//TODO 2 Matrix Transpose
//TODO 3 Hadamard Multiplication
//TODO 4 Matrix Slices / Views
//TODO 5 Element Wise Transformer Kernel

typedef enum{
    Double,
    Int,
    TypeCanary, /* <- for type checking asserts */
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
    UnsupportedTypeError,
    NoMatrixError,
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

#define let __auto_type


static inline MatrixResult matrix_create(MatrixType matrix_type,size_t rows,size_t cols){
    Matrix matrix;
    MatrixResult result;
    size_t matrix_size = 0;
    if(matrix_type == Double){
        matrix_size = sizeof(double);
    }
    else if(matrix_type == Int){
        matrix_size = sizeof(int64_t);
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

static inline void matrix_fill_int(Matrix matrix,int64_t value){
    int64_t* restrict matrix_data_int = (int64_t*)matrix.data;

    for_each(i,matrix){
        matrix_data_int[i] = value;
    }

}

static inline void matrix_fill_random_int(Matrix matrix,int64_t low,int64_t high){
    int64_t* restrict matrix_data_int = (int64_t*)matrix.data;

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



/* static const MatrixBinaryOpKernel matrix_add_kernels[] = {
    [Int]    = _matrix_add_int_kernel,
    [Double] = _matrix_add_double_kernel
};

static const MatrixBinaryOpKernel matrix_sub_kernels[] = {
    [Int]    = _matrix_sub_int_kernel,
    [Double] = _matrix_sub_double_kernel
}; */


static inline void _matrix_add_int_kernel(
    void* restrict _dst,
    const void* restrict _lhs,
    const void* restrict _rhs,
    size_t   elements_count
){
    int64_t* restrict dst = (int64_t*)_dst;
    int64_t* restrict lhs = (int64_t*)_lhs;
    int64_t* restrict rhs = (int64_t*)_rhs;

    for(size_t i = 0 ; i < elements_count ; ++i)
        dst[i] = lhs[i] + rhs[i];
}

static inline void _matrix_sub_int_kernel(
    void* restrict _dst,
    const void* restrict _lhs,
    const void* restrict _rhs,
    size_t   elements_count
){
    int64_t* restrict dst = (int64_t*)_dst;
    int64_t* restrict lhs = (int64_t*)_lhs;
    int64_t* restrict rhs = (int64_t*)_rhs;

    for(size_t i = 0 ; i < elements_count ; ++i)
        dst[i] = lhs[i] - rhs[i];
}

static inline void _matrix_mul_int_kernel(
    void* restrict _dst,
    const void* restrict _lhs,
    const void* restrict _rhs,
    size_t lhs_rows,
    size_t lhs_cols,
    size_t rhs_cols
){
    int64_t* restrict dst = (int64_t*)_dst;
    int64_t* restrict lhs = (int64_t*)_lhs;
    int64_t* restrict rhs = (int64_t*)_rhs;

    for(size_t i = 0 ; i < lhs_rows ; ++i){
        for(size_t j = 0 ; j < rhs_cols; ++j){
            int64_t sum = 0;
            for(size_t k = 0 ; k < lhs_cols; ++k){
                sum += lhs[i * lhs_cols + k] * rhs[k * rhs_cols + j];
            }
            dst[i * lhs_cols + j] = sum;
        }
    }
}

static inline void _matrix_clone_int_kernel(
    void* restrict _dst,
    const void* restrict _src,
    size_t rows,
    size_t cols 
){
    int64_t* restrict src_matrix_data = (int64_t*)_src;
    int64_t* restrict dst_matrix_data = (int64_t*)_dst;
        
    for(size_t i = 0 ; i < rows * cols ; ++i){
        dst_matrix_data[i] = src_matrix_data[i];
    }
}

static inline void _matrix_print_int_kernel(
    const void* restrict _lhs,
    size_t lhs_rows,
    size_t lhs_cols
){
    int64_t* lhs = (int64_t*)_lhs;
    for(size_t i = 0 ; i < lhs_rows; ++i){
        for(size_t j = 0 ; j < lhs_cols; ++j){
            printf("%3ld ",lhs[i * lhs_cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

static inline void _matrix_add_double_kernel(
    void* restrict _dst,
    const void* restrict _lhs,
    const void* restrict _rhs,
    size_t   elements_count
){
    double* restrict dst = (double*)_dst;
    double* restrict lhs = (double*)_lhs;
    double* restrict rhs = (double*)_rhs;
    for(size_t i = 0 ; i < elements_count ; ++i)
        dst[i] = lhs[i] + rhs[i];
}

static inline void _matrix_sub_double_kernel(
    void* restrict _dst,
    const void* restrict _lhs,
    const void* restrict _rhs,
    size_t   elements_count
){
    double* restrict dst = (double*)_dst;
    double* restrict lhs = (double*)_lhs;
    double* restrict rhs = (double*)_rhs;
    for(size_t i = 0 ; i < elements_count ; ++i)
        dst[i] = lhs[i] - rhs[i];
}

static inline void _matrix_mul_double_kernel(
    void* restrict _dst,
    const void* restrict _lhs,
    const void* restrict _rhs,
    size_t lhs_rows,
    size_t lhs_cols,
    size_t rhs_cols
){
    double* restrict dst = (double*)_dst;
    double* restrict lhs = (double*)_lhs;
    double* restrict rhs = (double*)_rhs;

    for(size_t i = 0 ; i < lhs_rows ; ++i){
        for(size_t j = 0 ; j < rhs_cols; ++j){
            double sum = 0;
            for(size_t k = 0 ; k < lhs_cols; ++k){
                sum += lhs[i * lhs_cols + k] * rhs[k * rhs_cols + j];
            }
            dst[i * lhs_cols + j] = sum;
        }
    }
}

static inline void _matrix_clone_double_kernel(
    void* restrict _dst,
    const void* restrict _src,
    size_t rows,
    size_t cols 
){
    double* restrict src_matrix_data = (double*)_src;
    double* restrict dst_matrix_data = (double*)_dst;
        
    for(size_t i = 0 ; i < rows * cols ; ++i){
        dst_matrix_data[i] = src_matrix_data[i];
    }
}

static inline void _matrix_print_double_kernel(
    const void* restrict _lhs,
    size_t lhs_rows,
    size_t lhs_cols
){
    double* lhs = (double*)_lhs;
    for(size_t i = 0 ; i < lhs_rows; ++i){
        for(size_t j = 0 ; j < lhs_cols; ++j){
            printf("%7.2f ",lhs[i * lhs_cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

typedef void(*MatrixBinaryOpKernel)( /* for addition and subtraction */
    void* restrict dst,
    const void* restrict lhs,
    const void* restrict rhs,
    size_t elements_count
);

typedef void(*MatrixUnaryOpGeneratingKernel)(
    void* restrict dst,
    const void* restrict lhs,
    size_t rows,
    size_t cols
);

typedef void(*MatrixUnaryOpNonGeneratingKernel)(
    const void* restrict lhs,
    size_t rows,
    size_t cols
);

typedef void(*MatrixMultiplicationKernel)( /* for multiplication as the name suggests*/
    void* restrict dst,
    const void* restrict lhs,
    const void* restrict rhs,
    size_t lhs_rows,
    size_t lhs_cols,
    size_t rhs_cols
);

typedef struct{
    MatrixBinaryOpKernel             add;
    MatrixBinaryOpKernel             sub;
    MatrixMultiplicationKernel       mul;
    MatrixUnaryOpGeneratingKernel    clone;
    MatrixUnaryOpNonGeneratingKernel print;
    //MatrixMultiplicationKernel div; // Multiplying with inverse
}MatrixVirtualTable;

static const MatrixVirtualTable matrix_ops_lookup[] ={
    [Int] = {
        .add    = _matrix_add_int_kernel,
        .sub    = _matrix_sub_int_kernel,
        .mul    = _matrix_mul_int_kernel,
        .clone  = _matrix_clone_int_kernel,
        .print  = _matrix_print_int_kernel,
    },
    [Double] = {
        .add    = _matrix_add_double_kernel,
        .sub    = _matrix_sub_double_kernel,
        .mul    = _matrix_mul_double_kernel,
        .clone  = _matrix_clone_double_kernel, 
        .print  = _matrix_print_double_kernel,
    }
};

static inline MatrixErrors
matrix_check_unary(const Matrix matrix){ // <- validates the matrix
    if(matrix.type >= TypeCanary){
        return UnsupportedTypeError;
    }
    return NoMatrixError;
}

static inline MatrixErrors 
matrix_check_binary(const Matrix lhs,const Matrix rhs) /* <- types,dimensions validator function */
{
    if(lhs.rows != rhs.rows || lhs.cols != rhs.cols){
        return InvalidDimensionError;
    }
    else if(lhs.type != rhs.type){
        return TypeMismatchError;
    }
    else if(lhs.type >= TypeCanary || rhs.type >= TypeCanary){
        return UnsupportedTypeError;
    }
    return NoMatrixError;
}

static inline MatrixErrors
matrix_check_multiplication(const Matrix lhs,const Matrix rhs)
{
    if(lhs.cols != rhs.rows){
        return InvalidDimensionError;
    }
    else if(lhs.type != rhs.type){
        return TypeMismatchError;
    }
    else if(lhs.type >= TypeCanary || rhs.type >= TypeCanary){
        return UnsupportedTypeError;
    }
    return NoMatrixError;
}

static inline MatrixResult
matrix_binary(const Matrix lhs,
              const Matrix rhs,
              MatrixBinaryOpKernel operation)
{
    MatrixResult result;
    let err = matrix_check_binary(lhs,rhs);

    if(err != NoMatrixError){   
        result.ok = false;
        result.matrix_error = err;
        return result;
    }

    result = matrix_create(lhs.type,lhs.rows,lhs.cols);
    if(result.ok == false)
        return result;

    size_t element_count = lhs.rows * lhs.cols;
    
    operation(result.matrix.data,
        lhs.data,
        rhs.data,
        element_count);
    return result;
}


static inline MatrixResult 
matrix_unary_generating(const Matrix matrix, /* <- an unary operation which takes self and nothing else which generates a new matrix,ex matrix_clone,matrix_transpose */
             MatrixUnaryOpGeneratingKernel operation) 
{
    MatrixResult result;
    let matrix_error = matrix_check_unary(matrix);

    if (matrix_error != NoMatrixError){
        result.ok = false;
        result.matrix_error = matrix_error;
        return result;
    }

    result = matrix_create(matrix.type,matrix.rows,matrix.cols);
    if(result.ok == false)
        return result;
    
    operation(result.matrix.data,
        matrix.data,
        matrix.rows,
        matrix.cols);
    return result;
}

static inline void 
matrix_unary_non_generating(const Matrix matrix, /* <- an unary operation which doesnt generate anything,like matrix_print */
             MatrixUnaryOpNonGeneratingKernel operation) 
{
    matrix_check_unary(matrix);

    operation(
        matrix.data,
        matrix.rows,
        matrix.cols);
}

static inline MatrixResult matrix_add(const Matrix lhs,const Matrix rhs){
    return matrix_binary(lhs,
        rhs,
        matrix_ops_lookup[lhs.type].add);
}

static inline MatrixResult matrix_sub(const Matrix lhs,const Matrix rhs){
    return matrix_binary(lhs,
        rhs,
        matrix_ops_lookup[lhs.type].sub);
}

static inline MatrixResult matrix_mul(const Matrix lhs,const Matrix rhs){
    MatrixResult result;

    let matrix_error = matrix_check_multiplication(lhs,rhs);
    if(matrix_error != NoMatrixError){
        result.ok = false;
        result.matrix_error = matrix_error;
        return result;    
    }
    
    result = matrix_create(lhs.type,lhs.rows,rhs.cols);
    if(result.ok == false)
        return result;

    matrix_ops_lookup[lhs.type].mul(result.matrix.data,
        lhs.data,
        rhs.data,
        lhs.rows,
        lhs.cols,rhs.cols);

    return result;
}

static inline void matrix_print(const Matrix matrix){
    matrix_unary_non_generating(matrix,
        matrix_ops_lookup[matrix.type].print);
}

static inline MatrixResult matrix_clone(const Matrix matrix){
   return matrix_unary_generating(matrix,
        matrix_ops_lookup[matrix.type].clone);
}

static inline bool matrix_equal(const Matrix lhs, const Matrix rhs){
    if(lhs.type != rhs.type)
        return false;
    if((lhs.rows != rhs.rows) || (lhs.cols != rhs.cols))
        return false;
    
    if(lhs.type == Int){
        int64_t* restrict lhs_data_int = (int64_t*)lhs.data;
        int64_t* restrict rhs_data_int = (int64_t*)rhs.data;
        
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

    drop Matrix m1 =  matrix_create(Int,3,3).matrix; // sort of like unwrap() <- use with extreme caution and never in production
    matrix_fill_int(m1,1);
    
    drop Matrix m2 = matrix_create(Double,2,3).matrix;
    matrix_fill_random_double(m2,1.0,9.0);


    drop Matrix m3 = matrix_clone(m1).matrix;
    drop Matrix m4 = matrix_clone(m2).matrix;

    return 0;
}
 


