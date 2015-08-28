

def YT(t,an,bn,trun=30, T=1,y0=5): # from 0 to t
    '''
    :param t: The point where the function need to be evaluated
    :param an: List/np.ndarray of coefficients of sin in Lambda, from a0 to an
    :param bn: List/np.ndarray of coefficients of cos in Lambda, from b0 to bn
    :param trun: Number of truncation terms
    :param T: Period
    :param y0: y(0)
    :return: function value
    '''
    trunAn=trunBn=trun

    if trunAn>len(an)-1:
        trunAn=len(an)-1

    if trunBn>len(bn)-1:
        trunBn=len(bn)-1
    an0=an[0]*0.5*t
    an_coe=np.asarray(an[1:])
    bn_coe=np.asarray(bn[1:])
    an_ran=np.arange(1,trunAn+1,dtype=float)
    bn_ran=np.arange(1,trunBn+1,dtype=float)
    an_val=an_coe/(2*np.pi*an_ran/T)*np.sin(2*np.pi*an_ran*t/T)
    bn_val=bn_coe/(2*np.pi*bn_ran/T)*(1-np.cos(2*np.pi*bn_ran*t/T))

    return y0*np.exp(an0+np.sum(an_val)+np.sum(bn_val))




def funOrigin(x,case=1):
    '''
    :param x:The point where the function need to be evaluated
    :param case: Considered case
    :return:Corresponding function value
    '''
    if case==1:
        return 32*x*x*(1-x)*(1-x)
    if case==2:
        return 20*x*(1-x)*(0.5-x)
    if case==3:
        return (1-x)*(x)
    if case==4:
        return x**3








def getcoef(expr, trun=30,T=1):
    '''
    :param expr: Symbolic expression of a function on which Fourier expansion is performed
    :param trun: Truncation terms
    :param T: Period
    :return: Fourier coefficients from 0 to trun
    '''
    an=([2/T*sympy.Integral(expr*sympy.cos(2*i*sympy.pi*x/T),(x,0,T)).evalf() for i in range(trun)])
    bn=([2/T*sympy.Integral(expr*sympy.sin(2*i*sympy.pi*x/T),(x,0,T)).evalf() for i in range(trun)])
    return an, bn


def coef(case=1):
    '''
    :param case: Considered case
    :return:Corresponding Fourier coefficients of trun=30
    '''
    if case==1:#32x^2(1-x)^2
        an=[2.133333333,-0.985534296,-0.061595894,-0.01216709,-0.003849743,-0.001576855,-0.000760443,-0.000410468,-0.000240609,-0.000150211,-9.86E-05,-6.73E-05,-4.75E-05,-3.45E-05,-2.57E-05,-1.95E-05,-1.50E-05,-1.18E-05,-9.39E-06,-7.56E-06,-6.16E-06,-5.07E-06,-4.21E-06,-3.52E-06,-2.97E-06,-2.52E-06,-2.16E-06,-1.85E-06,-1.60E-06,-1.39E-06,-1.22E-06]
        bn=np.zeros(31)
    if case==2:#20x(1-x)(0.5-x)
        an=np.zeros(31)
        bn=[0,0.967546033,0.120943254,0.035835038,0.015117907,0.007740368,0.00447938,0.002820834,0.001889738,0.001327224,9.68E-04,7.27E-04,5.60E-04,4.40E-04,3.53E-04,2.87E-04,2.36E-04,1.97E-04,1.66E-04,1.41E-04,1.21E-04,1.04E-04,9.09E-05,7.95E-05,7.00E-05,6.19E-05,5.50E-05,4.92E-05,4.41E-05,3.97E-05,3.58E-05]
    if case==3:#x(1-x)
        an=[0.333333333,-0.101321184,-0.025330296,-0.011257909,-0.006332574,-0.004052847,-0.002814477,-0.002067779,-0.001583143,-0.001250879,-0.001013212,-0.000837365,-0.000703619,-0.000599534,-0.000516945,-0.000450316,-0.000395786,-0.000350592,-0.00031272,-0.000280668,-0.000253303,-0.000229753,-0.000209341,-0.000191533,-0.000175905,-0.000162114,-0.000149883,-0.000138987,-0.000129236,-0.000120477,-0.000112579]
        bn=np.zeros(31)
    if case==4:#x^3
        an=[0.5,0.151981775,0.037995444,0.016886864,0.009498861,0.006079271,0.004221716,0.003101669,0.002374715,0.001876318,0.001519818,0.001256048,0.001055429,0.0008993,0.000775417,0.000675475,0.000593679,0.000525888,0.00046908,0.000421002,0.000379954,0.00034463,0.000314012,0.0002873,0.000263857,0.000243171,0.000224825,0.00020848,0.000193854,0.000180716,0.000168869]
        bn=[0,-0.269932585,-0.15310778,-0.104311543,-0.078821576,-0.063274959,-0.052827679,-0.045331799,-0.039694249,-0.035301404,-0.031782611,-0.028900916,-0.026497828,-0.024463356,-0.02271879,-0.021206325,-0.019882557,-0.018714264,-0.017675587,-0.016746099,-0.015909447,-0.01515239,-0.014464088,-0.013835584,-0.013259412,-0.012729299,-0.012239935,-0.011786797,-0.011366006,-0.010974219,-0.010608538]
    return np.asarray(an),np.asarray(bn)


def YTOrigin(t,y0=5,case=1):
    '''
    :param t: The point where the function need to be evaluated
    :param y0: y(0)
    :param case: Considered case
    :return:Function value
    '''
    return y0*np.exp(integrate.quad(funOrigin, 0, t,args=(case))[0])





def coef_random(an,bn,Random_trun=30,E=1,p=-1):
    '''
    takes Fourier coefficients and return a set of random coefficients in a cube around the input.
    :param an, bn:Base Fourier coefficient
    :param Random_trun:Until which term the random perturbation will be added
    :param E:Half width of interval
    :param p:Decay rate
    :return:Random coefficients an, bn
    '''
    bn_ran=an_ran=np.arange(1,Random_trun+2,dtype=float)
    an_random=an_ran**p*E*np.random.uniform(-1,1,(Random_trun+1,))
    bn_random=bn_ran**p*E*np.random.uniform(-1,1,(Random_trun+1,))
    an_random.resize(an.shape)
    bn_random.resize(bn.shape)
    an_random=an+an_random
    bn_random=bn+bn_random
    return an_random,bn_random


def getint_random(an,bn,Random_trun=30,E=1,p=-1,a=0,b=1):
    '''
    Takes Fourier coefficients and output random integration value
    :param an, bn: Base Fourier coefficients
    :param Random_trun:until which term the perturbation will be added
    :param E:Half width of interval
    :param p:Decay rate
    :param a:Lower limit of integration
    :param b:Upper limit of integration
    :return:Integration value
    '''
    an_random,bn_random=coef_random(an,bn,Random_trun=Random_trun,E=E,p=p)
    int_random=integrate.quad(YT, a, b,args=(an_random,bn_random))[0]
    return int_random
