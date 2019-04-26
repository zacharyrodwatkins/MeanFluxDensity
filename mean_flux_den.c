//Gotta make the flux
#include <python2.7/Python.h>
#include <pthread.h>
#include "photometry.h"
pthread_mutex_t lock;


static void update_sum(double *sum , double operand){
    pthread_mutex_lock(&lock); //Obtain lock to update sum
    *sum += operand;
    pthread_mutex_unlock(&lock);
}


static void get_doubles(PyObject* lambda_py , PyObject* trans_py, PyObject* flux_py,
double* lambda , double* trans, double* flux, int index){
        
        PyObject *lamb_pyf, *flux_pyf, *trans_pyf;
        PyObject *lambda_i = PySequence_Fast_GET_ITEM(lambda_py, index);
        PyObject *flux_lamb_i = PySequence_Fast_GET_ITEM(flux_py, index);
        PyObject *trans_i = PySequence_Fast_GET_ITEM(trans_py, index);

        lamb_pyf  = PyNumber_Float(lambda_i);
        flux_pyf = PyNumber_Float(flux_lamb_i);
        trans_pyf = PyNumber_Float(trans_i);
        
        if (!lamb_pyf||!flux_pyf||!trans_pyf) {
            Py_DECREF(lambda_py);
            Py_DECREF(flux_py);
            Py_DECREF(trans_py);
            PyErr_SetString(PyExc_TypeError, "all items must be numbers");
        }
        
        *lambda = PyFloat_AS_DOUBLE(lamb_pyf);
        *flux = PyFloat_AS_DOUBLE(flux_pyf);
        *trans = PyFloat_AS_DOUBLE(trans_pyf);

        Py_DECREF(lamb_pyf);
        Py_DECREF(flux_pyf);
        Py_DECREF(trans_pyf);
}

//converts to doubles, computes integrals
static void sequent_work(int start, int stop, PyObject* lambda_py, PyObject* Trans_py,
 PyObject* flux_lambda_py, double* total_sum , double* norm_sum){
 
        double *prev_lamb = malloc(sizeof(double));
        double *prev_flux = malloc(sizeof(double));
        double *prev_trans_p = malloc(sizeof(double));
        double *i_lamb = malloc(sizeof(double));
        double *i_flux = malloc(sizeof(double));
        double *i_trans = malloc(sizeof(double));

        get_doubles(lambda_py, Trans_py, flux_lambda_py,
            prev_lamb, prev_trans_p, prev_flux, start);

        double flux_sum_loc = 0;
        double norm_sum_loc = 0;
        double f_nu_prev = get_f_nu(*prev_flux*1e-8 /*1e-8 factor converts from ergs/s/cm^/cm to ergs/s/cm^/A*/, *prev_lamb);
        double nu_prev = get_nu(*prev_lamb);
        double f_nu;
        double nu;
        double prev_trans = *prev_trans_p;
        free(prev_trans_p);
        free(prev_lamb);
        free(prev_flux);
        for (int i = start+1; i < stop; i++){
            get_doubles(lambda_py, Trans_py, flux_lambda_py,
            i_lamb, i_trans, i_flux, i);
            f_nu = get_f_nu(*i_flux*1e-8 /*1e-8 factor converts from ergs/s/cm^/cm to ergs/s/cm^/A*/, *i_lamb);
            nu = get_nu(*i_lamb);
            flux_sum_loc += (nu_prev-nu)*0.5*
                (flux_den(f_nu, *i_trans, nu)+flux_den(f_nu_prev, prev_trans, nu_prev));
            norm_sum_loc += (nu_prev-nu)*0.5*
                (filter_norm(*i_trans, nu)+filter_norm(prev_trans, nu_prev));

            prev_trans = *i_trans;
            f_nu_prev = f_nu;
            nu_prev = nu;
        }        
        update_sum(total_sum, flux_sum_loc);
        update_sum(norm_sum, norm_sum_loc*C);
        free(i_flux);
        free(i_trans);
        free(i_lamb); 
 }


/* here is how you expose it to Python code: */
static PyObject *totalDoubles(PyObject *self, PyObject *args)
{   

    pthread_mutex_init(&lock, NULL);
    PyObject* lambda, * trans,* flux_lambda;
    double *flux_sum = malloc(sizeof(double));
    *flux_sum = 0;
    double *norm_sum = malloc(sizeof(double));
    *norm_sum = 0;
    double result;
    int len;


    /* get one argument as a sequence */
    if(!PyArg_ParseTuple(args, "OOO", &lambda, &trans, &flux_lambda))
        return 0;

    flux_lambda = PySequence_Fast(flux_lambda, "argument must be iterable");
    trans = PySequence_Fast(trans, "argument must be iterable");
    lambda = PySequence_Fast(lambda, "argument must be iterable");
    if(!flux_lambda||!trans||!lambda)
        return 0;

    /* prepare data as an array of doubles */
    len = PySequence_Fast_GET_SIZE(lambda);

    /* converts the lambda trans and flux_lambda python objects to doubles and computes 
    *  the intergrals
    * */
    sequent_work(0 , len, lambda, trans, flux_lambda, flux_sum, norm_sum);
   
    /* clean up, compute, and return result */
    Py_DECREF(lambda);
    Py_DECREF(flux_lambda);
    Py_DECREF(trans);
    result = (*flux_sum)/(*norm_sum);
    //free(flux_lambda);
    //free(trans);
    //free(lambda);
    free(flux_sum);
    free(norm_sum);
    return Py_BuildValue("d", result);
}

static PyMethodDef totalMethods[] = {
    {"total", totalDoubles, METH_VARARGS, "Sum a sequence of numbers."},
    {0} /* sentinel */
};

void
initFluxDen (void) {
     (void) Py_InitModule("FluxDen", totalMethods);
}


int main(void) {
    return 0;
}


