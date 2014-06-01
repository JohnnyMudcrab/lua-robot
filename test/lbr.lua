local robot = require 'robot'
local link = require 'link'
local matrix = require 'matrix'

pi12 = math.pi * 0.5

a     = {0, 0, 0, 0, 0, 0, 0}
alpha = {0, pi12, -pi12, pi12, -pi12, pi12, -pi12}
d     = {0.2, 0, 0.4, 0, 0.39, 0, 0.0675}

m1 = 2.7;
m2 = 2.7;
m3 = 2.7;
m4 = 2.7;
m5 = 1.7;
m6 = 1.6;
m7 = 0.3;

rc1 = matrix {0.00015950202941895, -0.025668799877167, -0.080369353294373};
rc2 = matrix {-3.4362077713013e-05, 0.0888432264328, -0.027725577354431};
rc3 = matrix {-5.4091215133667e-05, 0.027707517147064, -0.08038717508316};
rc4 = matrix {2.7269124984741e-05, 0.08891624212265, 0.025711119174957};
rc5 = matrix {3.361701965332e-05, -0.025357067584991, -0.095635056495667};
rc6 = matrix {0.00021129846572876, -0.0057910680770874, -0.0088013410568237};
rc7 = matrix {-1.189112663269e-05, -0.00049149990081787, 2.1457672119141e-05};

I1 = matrix {{0.016340516507626,	0,	0},
             {0,	0.016173135489225,	0.003534},
             {0,	0.003534,	0.0050259656272829}}

I2 = matrix {{0.01634038425982,	0,	0},
             {0,	0.016173111274838,	0.003534},
             {0,	0.003534,	0.0050272806547582}}

I3 = matrix {{0.016340393573046,	0,	0},
             {0,	0.016172878444195,	-0.003533},
             {0,	-0.003533,	0.0050275097601116}}

I4 = matrix {{0.01634038053453,	0,	0},
             {0,	0.016172980889678,	-0.003534},
             {0,	-0.003534,	0.0050261416472495}}

I5 = matrix {{0.009819190017879,	0,	0},
             {0,	0.0090915467590094,	0.003093},
             {0,	0.003093,	0.0037081523332745}}

I6 = matrix {{0.0030111982487142,	0,	0},
             {0,	0.003022399963811,	0},
             {0,	0,	0.0034143980592489}}

I7 = matrix {{0.0001017179529299,	0,	0},
             {0,	0.00010171789472224,	0},
             {0,	0,	0.00015843522851355}}


L = {}
L[1] = link.new(a[1], alpha[1], d[1], m1, rc1, I1)
L[2] = link.new(a[2], alpha[2], d[2], m2, rc2, I2)
L[3] = link.new(a[3], alpha[3], d[3], m3, rc3, I3)
L[4] = link.new(a[4], alpha[4], d[4], m4, rc4, I4)
L[5] = link.new(a[5], alpha[5], d[5], m5, rc5, I5)
L[6] = link.new(a[6], alpha[6], d[6], m6, rc6, I6)
L[7] = link.new(a[7], alpha[7], d[7], m7, rc7, I7)

lbr = robot.new(L)

print('\nDH-Parameter')
lbr:printDH()

q1 = {0, 0, 0, 0, 0, 0, 0};
q2 = {0, 0.3491, 0, 0.3491, 0, 0, 0};

qd = {1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3}

qdd = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5}

grav = matrix {{0},{0},{-9.81}}
tau = lbr:gravLoad(q2,grav)

print('\ngravLoad')
for i = 1,7,1 do
    print(tau[i])
end

C = lbr:coriolis(q2,qd)

print('\ncoriolis')
matrix.print(C)
