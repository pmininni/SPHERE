! Initial magnetic field
         q1 = 3
         l1 = 3
         xsib(l1*(l1+1)/2,q1) = b0*1.d0
         xsib(l1*(l1+1)/2+1,q1) = b0*(1.d0-IM)
         xsib(l1*(l1+1)/2+2,q1) = b0*(1.d0-IM)
         xsib(l1*(l1+1)/2,2*q-q1+1) = .6*b0*1.d0
         xsib(l1*(l1+1)/2+1,2*q-q1+1) = .6*b0*(1.d0-IM)
         xsib(l1*(l1+1)/2+2,2*q-q1+1) = .6*b0*(1.d0-IM)
