p=419;K.<a>=QuadraticField(-p);E=EllipticCurve(GF(419), [0,0,0,1,0])
hk=K.class_number()
r=E.order();print(r)
def l_Point(l,E0,p):
    cofactor=(p+1)/l
    while True:
        Q = E0.random_point()
        P = cofactor * Q                
        if P.order() == l:
            break
    return P
def walk(l,E,p):
    P = E(0).division_points(l)[1]             
    return E.isogeny(P).codomain()
def Back(l,e,A,p):
    F=GF(p);d=F(2)
    Et=EllipticCurve(F, [0,A*d,0,d^2,0])
    if e<0:
        e=(-e)%hk
    for _ in range(abs(e)):
        P=l_Point(l,Et,p)
        Et=Et.isogeny(P).codomain()
        M=Et.montgomery_model();u2 = int(sqrt(F(1)/(M.a4())))
    return [(M.a2())*u2, -((M.a2())*u2)]
for A in Back(3,-1,0,p):
    EA = EllipticCurve(GF(p), [0, A, 0, 1, 0])
    if walk(3,EA,p).is_isomorphic(E):
        print(A)
