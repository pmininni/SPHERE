! Initial velocity field
         q1 = 3
         l1 = 3
         xsiv(l1*(l1+1)/2,q1) = -u0*1.d0
         xsiv(l1*(l1+1)/2+1,q1) = u0*(1.d0+IM)
         xsiv(l1*(l1+1)/2+2,q1) = u0*(1.d0+IM)
         xsiv(l1*(l1+1)/2,2*q-q1+1) = -u0*1.d0
         xsiv(l1*(l1+1)/2+1,2*q-q1+1) = u0*(1.d0+IM)
         xsiv(l1*(l1+1)/2+2,2*q-q1+1) = u0*(1.d0+IM)
