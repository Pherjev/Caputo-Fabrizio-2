import numpy
import matplotlib.pyplot as plt

L = 0.01
gamma = 0.01
gamma1 = 0.01
gamma2 = 0.05
k = 0.1
mu = 0.001
alpha = 0.12
delta = 0.001
q = 2


# Parametros fraccionales

def M(x):
        return 2./(2.-x)

sigma = 0.99

M1 = (2*(1-sigma))/float(M(sigma)*(2-sigma))

M2 = (2*sigma)/float((2-sigma)*M(sigma))



# Funciones simples

def f(S,I):
        return k*I*S**q

def g1(S,I):
	return k*q*I*(S**(q-1))

def g2(S,I):
	return k*(S**q)


def functionS(S,I,R):
        return L - mu*S - f(S,I) + gamma2*I + delta*R

def functionI(S,I,R):
        return f(S,I) -(mu + gamma1 + gamma2 + alpha)*I

def functionR(S,I,R):
        return gamma2*I - (mu + delta)*R


# Funciones fraccionarias

def c1(S,I,R):
	return -mu-g1(S,I)

def c2(S,I,R):
	return -g2(S,I) + gamma1

def c3(S,I,R):
	return g1(S,I)

def c4(S,I,R):
	return g2(S,I) - (mu + gamma1 + gamma2 + alpha)

c5 =  -mu - delta


def k1(S,I,R):
	return (M1*c2(S,I,R))/float(1-M1*c1(S,I,R))

def k2(S,I,R):
        return (M1*delta)/float(1-M1*c1(S,I,R))

def k3(S,I,R):
        return (M2*functionS(S,I,R))/float(1-M1*c1(S,I,R))

def k4(S,I,R):
        return (M1*c3(S,I,R))/float(1-M1*c4(S,I,R))

def k5(S,I,R):
        return (M2*functionI(S,I,R))/float(1-M1*c4(S,I,R))

def k6(S,I,R):
        return (M1*gamma2)/float(1-M1*c5)

def k7(S,I,R):
        return (M2*functionR(S,I,R))/float(1-M1*c5)



def dfunctionS(S,I,R):
	num = k1(S,I,R)*k5(S,I,R) + k2(S,I,R)*k5(S,I,R)*k6(S,I,R) + k2(S,I,R)*k7(S,I,R) + k3(S,I,R)
	den = 1 - k1(S,I,R)*k4(S,I,R) - k2(S,I,R)*k4(S,I,R)*k6(S,I,R)
	return num/float(den) 

def dfunctionI(S,I,R):
	return k4(S,I,R)*dfunctionS(S,I,R) + k5(S,I,R)

def dfunctionR(S,I,R):
	return k6(S,I,R)*dfunctionI(S,I,R) + k7(S,I,R)

h = 0.001

s = 0.9
i = 0.09
r = 0.01
t = 0

St = [s]
It = [i]
Rt = [r]
T  = [t]

for j in range(400000):

	t += h
	s1 = s + h*functionS(s,i,r)
	i1 = i + h*functionI(s,i,r)
	r1 = r + h*functionR(s,i,r)

	s = s1 + 0
	i = i1 + 0
	r = r1 + 0

	St.append(s)
	It.append(i)
	Rt.append(r)
	T.append(t)

sigma = 0.9

M1 = (2*(1-sigma))/float(M(sigma)*(2-sigma))

M2 = (2*sigma)/float((2-sigma)*M(sigma))



s = 0.9
i = 0.09
r = 0.01
t = 0

St2 = [s]
It2 = [i]
Rt2 = [r]
T2  = [t]

for j in range(400000):

        t += h
        s1 = s + h*dfunctionS(s,i,r)
        i1 = i + h*dfunctionI(s,i,r)
        r1 = r + h*dfunctionR(s,i,r)

        s = s1 + 0
        i = i1 + 0
        r = r1 + 0

        St2.append(s)
        It2.append(i)
        Rt2.append(r)
        T2.append(t)

sigma = 0.6

M1 = (2*(1-sigma))/float(M(sigma)*(2-sigma))

M2 = (2*sigma)/float((2-sigma)*M(sigma))



s = 0.9
i = 0.09
r = 0.01
t = 0

St3 = [s]
It3 = [i]
Rt3 = [r]
T3  = [t]

for j in range(400000):

        t += h
        s1 = s + h*dfunctionS(s,i,r)
        i1 = i + h*dfunctionI(s,i,r)
        r1 = r + h*dfunctionR(s,i,r)

        s = s1 + 0
        i = i1 + 0
        r = r1 + 0

        St3.append(s)
        It3.append(i)
        Rt3.append(r)
        T3.append(t)


sigma = 0.3

M1 = (2*(1-sigma))/float(M(sigma)*(2-sigma))

M2 = (2*sigma)/float((2-sigma)*M(sigma))


s = 0.9
i = 0.09
r = 0.01
t = 0

St4 = [s]
It4 = [i]
Rt4 = [r]
T4  = [t]

for j in range(400000):

        t += h
        s1 = s + h*dfunctionS(s,i,r)
        i1 = i + h*dfunctionI(s,i,r)
        r1 = r + h*dfunctionR(s,i,r)

        s = s1 + 0
        i = i1 + 0
        r = r1 + 0

        St4.append(s)
        It4.append(i)
        Rt4.append(r)
        T4.append(t)


#print(St)

fig,ax = plt.subplots()

ax.plot(T,St, color = (0.4,0.7,0.3),label='1')
ax.plot(T,St2,ls='--', color = (0.4,0.7,0.3),label='0.9')
ax.plot(T,St3,ls='-.', color = (0.4,0.7,0.3),label='0.6')
ax.plot(T,St4,ls=':', color = (0.4,0.7,0.3), label='0.3')
#ax.plot(T,It)
#ax.plot(T,Rt)
#ax.set_title(u"SIRS Fractional Model")
ax.spines['left'].set_position(('outward',10))
ax.spines['bottom'].set_position(('outward',10))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_yticks(range(10),minor=True)
ax.set_ylim([0,2])
ax.legend(framealpha=1)
#ax.legend(bbox_to_anchor=(1.1, 1.05))
ax.legend()#loc="upper left", bbox_to_anchor=(0.8,0.2))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
fig.savefig('S1-2.png')
plt.show()


fig,ax = plt.subplots()

#print(It)

ax.plot(T,It, color = (0.7,0.3,0.3),label='1')
ax.plot(T,It2,ls='--', color = (0.7,0.3,0.3),label='0.9')
ax.plot(T,It3,ls='-.', color = (0.7,0.3,0.3),label='0.6')
ax.plot(T,It4,ls=':', color = (0.7,0.3,0.3), label='0.3')
#ax.plot(T,It)
#ax.plot(T,Rt)
#ax.set_title(u"SIRS Fractional Model")
ax.spines['left'].set_position(('outward',10))
ax.spines['bottom'].set_position(('outward',10))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_yticks(range(10),minor=True)
ax.set_ylim([0,2])
ax.legend(framealpha=1)
#ax.legend(bbox_to_anchor=(1.1, 1.05))
ax.legend()#loc="upper left", bbox_to_anchor=(0.8,0.2))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
fig.savefig('I1-2.png')
plt.show()

fig,ax = plt.subplots()

ax.plot(T,Rt, color = (0.3,0.4,0.7),label='1')
ax.plot(T,Rt2,ls='--', color = (0.3,0.4,0.7),label='0.9')
ax.plot(T,Rt3,ls='-.', color = (0.3,0.4,0.7),label='0.6')
ax.plot(T,Rt4,ls=':', color = (0.3,0.4,0.7), label='0.3')
#ax.plot(T,It)
#ax.plot(T,Rt)
#ax.set_title(u"SIRS Fractional Model")
ax.spines['left'].set_position(('outward',10))
ax.spines['bottom'].set_position(('outward',10))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_yticks(range(10),minor=True)
ax.set_ylim([0,2])
ax.legend(framealpha=1)
#ax.legend(bbox_to_anchor=(1.1, 1.05))
ax.legend()#loc="upper left", bbox_to_anchor=(0.8,0.2))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
fig.savefig('R1-2.png')
plt.show()
