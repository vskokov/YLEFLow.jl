module YLEFlow

using Symbolics

# Write your package code here.

@variables t k ρ a0 a1 a2 a3 a4 a5 a6 a7 a8
@variables U(..) af1(..) af2(..) af3(..) af4(..) af5(..) af6(..) af7(..) af8(..)

msigma2 = 2*ρ*Differential(ρ)(Differential(ρ)(U(ρ))) + Differential(ρ)(U(ρ))
#dUdk = exp(5*t)/(exp(2*t) + msigma2)
dUdk = k^5/(k^2 + msigma2)


S1 =  [U(ρ),af1(ρ),af2(ρ),af3(ρ),af4(ρ),af5(ρ),af6(ρ),af7(ρ),af8(ρ)]
S2 =  [af1(ρ),af2(ρ),af3(ρ),af4(ρ),af5(ρ),af6(ρ),af7(ρ),af8(ρ),0]
S3 = [a0, a1, a2, a3, a4, a5, a6, a7, a8]

R=[ρ,ρ,ρ,ρ,ρ,ρ,ρ,ρ,ρ]
R[1]=dUdk

Q=similar(R)
for (c,p) in zip(S1,S2)
    R[1]=simplify(substitute(R[1],Dict( Differential(ρ)(c)  => p)))
end

for i in 2:length(R)
    R[i] = Symbolics.derivative(R[i-1],ρ)

    for (c,p) in zip(S1,S2)
        R[i]=(substitute(R[i],Dict( Differential(ρ)(c)  => p)))
    end
end

Q = copy(R)

for (c,p) in zip(S1,S3)
    for i in 1:length(R)
        Q[i]=substitute(Q[i],Dict( c  => p))
    end
end

for i in 1:length(R)
    Q[i]=simplify(Q[i])
end

Symbolics.derivative(Q[2],a2)

dU1dk = Symbolics.build_function(Q[2],k,a1,a2,a3,a4,a5,a6,a7,a8,ρ,expression=Val{false})
dU2dk = Symbolics.build_function(Q[3],k,a1,a2,a3,a4,a5,a6,a7,a8,ρ,expression=Val{false})
dU3dk = Symbolics.build_function(Q[4],k,a1,a2,a3,a4,a5,a6,a7,a8,ρ,expression=Val{false})
dU4dk = Symbolics.build_function(Q[5],k,a1,a2,a3,a4,a5,a6,a7,a8,ρ,expression=Val{false})
dU5dk = Symbolics.build_function(Q[6],k,a1,a2,a3,a4,a5,a6,a7,a8,ρ,expression=Val{false})
dU6dk = Symbolics.build_function(Q[7],k,a1,a2,a3,a4,a5,a6,a7,a8,ρ,expression=Val{false})
dU7dk = Symbolics.build_function(Q[8],k,a1,a2,a3,a4,a5,a6,a7,a8,ρ,expression=Val{false})
dU8dk = Symbolics.build_function(Q[9],k,a1,a2,a3,a4,a5,a6,a7,a8,ρ,expression=Val{false})


function drodk(dU1dk,dU2dk,a2,a3,ρ)
    return  -(dU1dk + 2*ρ*dU2dk)/(3*a2+2*ρ*a3)
end

function flow!(da,a,m2,t)
 k = exp(t)
 ρ = a[1]
 a1 = m2 - 2ρ * a[2]

 dU1dk_s = dU1dk(t, a1, a[2], a[3], a[4], a[5], a[6], a[7], a[8], ρ)
 dU2dk_s = dU2dk(t, a1, a[2], a[3], a[4], a[5], a[6], a[7], a[8], ρ)
 dU3dk_s = dU3dk(t, a1, a[2], a[3], a[4], a[5], a[6], a[7], a[8], ρ)
 dU4dk_s = dU4dk(t, a1, a[2], a[3], a[4], a[5], a[6], a[7], a[8], ρ)
 dU5dk_s = dU5dk(t, a1, a[2], a[3], a[4], a[5], a[6], a[7], a[8], ρ)
 dU6dk_s = dU6dk(t, a1, a[2], a[3], a[4], a[5], a[6], a[7], a[8], ρ)
 dU7dk_s = dU7dk(t, a1, a[2], a[3], a[4], a[5], a[6], a[7], a[8], ρ)
 dU8dk_s = dU8dk(t, a1, a[2], a[3], a[4], a[5], a[6], a[7], a[8], ρ)


 da[1] = drodk(dU1dk_s,dU2dk_s,a[2],a[3],ρ)

 da[2] = dU2dk_s + a[3] * da[1]
 da[3] = dU3dk_s + a[4] * da[1]
 da[4] = dU4dk_s + a[5] * da[1]
 da[5] = dU5dk_s + a[6] * da[1]
 da[6] = dU6dk_s + a[7] * da[1]
 da[7] = dU7dk_s + a[8] * da[1]
 da[8] = dU8dk_s
end

function g()

end

end
