import matplotlib.pyplot as plt
import math as ma
import numpy as np
##################################### Диапазоны, шаги
diap_x = 1.
diap_y = 2.
diap_t = 2.
nx = 30.
ny = 30.
nt = 30.
############################################## система ввода значений
print("В каком временном интервале решать задачу (0,t): t =")
diap_t = float(input())
print("Количество шагов по x")
nx = int(input())
print("Количество шагов по y")
ny = int(input())
print("Количество шагов по t")
nt = int(input())
################################################## инициализация сетки и шагов
setka = np.zeros((nx,ny,2*nt+1), dtype= float) #сетка решений, 2nt + 1 потому что приходится работать с t+1/2/t-1/2
x_step = diap_x/float(nx-1)
y_step = diap_y/float(ny-1) #размеры шагов
t_step = diap_t/float(nt-1)
gamma_x = t_step/x_step/x_step
gamma_y = t_step/y_step/y_step
################################################начальные условия
def f_1(i,j,k): # функция для первого прогона.
    return -(0.5*gamma_y*(setka[i,j-1,k-1]+setka[i,j+1,k-1])+(1-gamma_y)*setka[i,j,k-1] +0.5*t_step*(t_step/2*(k+1)*(t_step/2*((k+1)))*(j*y_step)))
def f_2(i,j,k):# функция второго прогона
    return -(0.5*gamma_x*(setka[i-1,j,k-1]+setka[i+1,j,k-1])+(1-gamma_x)*setka[i,j,k-1]+0.5*t_step*(t_step/2*(k-1))*(t_step/2*(k-1)*(j*y_step)))
############################################################################
def progonka_1(j,k):             #функция прогона 1. По иксам
    alpha = np.zeros(nx)
    beta = np.zeros(nx)
    alpha[1] = 1
    beta[1] = 0
    A_x = 0.5*gamma_x
    B_x = 1+gamma_x
    C_x = 0.5* gamma_x
    for z in range(1, nx-1,1):
        F_x = f_1(z,j,k)
        alpha[z+1] = C_x/(B_x-A_x*alpha[z])
        beta [z+1]  = (-A_x*beta[z]+F_x)/(-B_x+A_x*alpha[z])
    setka[nx-1,j,k] = 0
    for z in range(nx-1,0, -1):
        setka[z-1,j,k] = alpha[z]*setka[z,j,k]+beta[z]
    setka[0,j,k] = setka[1,j,k]
###############################################################
def progonka_2(i,k):
    alpha = np.zeros(ny)
    beta = np.zeros(ny)
    alpha[1] = 1
    beta[1] = 0
    A_y = 0.5 * gamma_y
    B_y = 1+gamma_y
    C_y = 0.5* gamma_y
    for z in range(1, ny - 1, 1):
        F_x = f_2(i, z, k)
        alpha[z + 1] = C_y / (B_y - A_y * alpha[z])
        beta[z + 1] = (-A_y * beta[z] + F_x) / (-B_y + A_y * alpha[z])
    setka[i, ny-1, k] = (beta[-1])/(1 - alpha[-1])
    for z in range(ny - 1, 0, -1):
        setka[i, z-1, k] = alpha[z] * setka[i, z, k] + beta[z]
    setka[i,0,k] = setka[i,1,k]
#######################################################################
def plotting(t):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    b = np.zeros((nx, ny), dtype=float)
    for i in range(0, nx, 1):
        for j in range(0, ny, 1):
            b[i][j] = setka[i][j][int(t/diap_t*2*nt)]
    b = b.transpose()
    x = np.arange(0,diap_x+x_step/2,x_step)
    y = np.arange(0,diap_y+y_step/2,y_step)
    X, Y = np.meshgrid(x,y)
    #b = np.zeros((nx,ny),dtype = float)
    ax.plot_surface(X, Y, b, rstride=1, cstride=1,
                cmap='coolwarm', edgecolor='none')
    plt.show()
###########################################################
def graph_builder_y(t,x):
    fig = plt.figure()
    b = np.zeros(ny, dtype=float)
    for j in range(0, ny, 1):
        b[j] = setka[int(x/x_step),j,int(t/diap_t*2*nt)]
    b = b.transpose()
    y = np.arange(0,diap_y+y_step/2,y_step)
    plt.plot(y, b)
    plt.show()
#############################################################################
def graph_builder_x(t,y):
    fig = plt.figure()
    b = np.zeros(nx, dtype=float)
    for j in range(0, nx, 1):
        b[j] = setka[j,int(y/y_step),int(t/diap_t*2*nt)]
    b = b.transpose()
   # y = np.arange(0,diap_y+y_step/2,y_step)
    x = np.arange(0, diap_x + x_step / 2, x_step)
    plt.plot(x, b)
    plt.show()
################################################################################
for x in range(0,nx):
    for y in range(0,ny):
        setka[x,y,0] =(x*x_step*x_step*x-1)*ma.cos(ma.pi*y*y_step)
#############################################################################
for k in range(1,2*nt,2):
    for j in range(1,ny-1):
        progonka_1(j,k)
    for i in range(1,nx-1):
        progonka_2(i,k+1)
    setka[nx - 1, :, k + 1] = 0
    setka[0, :, k + 1] = setka[1, :, k + 1]
    setka[:, ny - 1, k + 1] = setka[:, ny - 2, k + 1]
    setka[:, 0, k ] = setka[:, 1, k ]
    setka[nx - 1, :, k ] = 0
    setka[0, :, k ] = setka[1, :, k ]
    setka[:, ny - 1, k ] = setka[:, ny - 2, k ]
    setka[:, 0, k ] = setka[:, 1, k ]


#graph_builder_y(diap_t,0.25)
#graph_builder_x(diap_t,0.25)
def build_surf_y(x):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_xlabel('y')
    ax.set_ylabel('t')
    ax.set_zlabel('u')
    b = np.zeros((nx, 2*nt+1), dtype=float)
    for j in range(0, ny, 1):
        for k in range(0, 2 * nt + 1, 1):
            b[j][k] = setka[int(x / diap_x *(nx - 1))][j][k]
    b = b.transpose()
    x = np.arange(0, diap_y + y_step / 2, y_step)
    t = np.arange(-t_step/2, diap_t + t_step / 2 +t_step/180, t_step / 2)
    Y, T = np.meshgrid(x, t)
    # b = np.zeros((nx,ny),dtype = float)
    ax.plot_surface(Y, T, b, rstride=1, cstride=1,
                    cmap='coolwarm', edgecolor='none')
    plt.show()

def build_surf_x(y):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u')
    b = np.zeros((nx, 2*nt+1), dtype=float)
    for i in range(0, nx, 1):
        for k in range(0, 2*nt+1, 1):
            b[i][k] = setka[i][int(y / diap_y*(ny-1))][k]
    b = b.transpose()
    x = np.arange(0, diap_x + x_step / 2, x_step)
    t = np.arange(-t_step/2, diap_t + t_step / 2+t_step/180, t_step / 2)
    print(x.size)
    print(t.size)
    print(b.size)
    X, T = np.meshgrid(x, t)
    # b = np.zeros((nx,ny),dtype = float)
    ax.plot_surface(X,T , b, rstride=1, cstride=1,
                    cmap='coolwarm', edgecolor='none')
    plt.show()


def theory_solution(n,m,x,y,t):
    u_sum = 0.0
    for i in range(0,n,1):
        for j in range(1,m,1):
            u = T(t,i,j)*ma.cos(ma.pi*(2*j-1)/2*x)*ma.cos(ma.pi*i*y/2)
            u_sum = u_sum+u

    return u_sum
def T(t,n,m):
    res = 0.0
    lam = float((ma.pi/2*(2*m-1))**2 + (ma.pi*n/2)**2)
    #print(lam)
    if(n == 2):
        res = 32*((-1)**m)/(ma.pi**3)/(2*n-1)**3*ma.exp(-lam*t) + 1/lam**3*(A(n,m))*(((lam*t)**2-2*lam*t+2)-2* ma.exp(-lam*t))
    else:
        res = 1/lam**3*(A(n,m))*(((lam*t)**2-2*lam*t+2)-2* ma.exp(-lam*t))
    return res
def A(n,m):
    a = 0.0
    if(n==0):
        a = 4 * (-1)**m/(ma.pi*(1 - 2*m))
    else:
        a = 16 *(-1)**m/(1-2*m)*((-1)**n-1)/(ma.pi**3*m**2)
    #print(a)
    return a

print("Введите x")
x = 0.25 #float(input())
print("Введите y")
y = 0.25 #float(input())
print("Введите t")
t = diap_t

ind_x = int(x / diap_x*(nx-1))
ind_y = int(y / diap_y*(ny-1))
ind_t = int(t/diap_t*2*nt)
#plotting(diap_t)
u = theory_solution(100,100,x,y,t)
#print(u)

print(setka[ind_x,ind_y,ind_t])
print(abs(u - setka[ind_x,ind_y,ind_t]))

#build_surf_x(1.5)
#build_surf_y(0.75)
plotting(diap_t)