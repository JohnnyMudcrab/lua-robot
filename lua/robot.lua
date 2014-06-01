-------------------------------------------------
--
-- robot.lua
--
-------------------------------------------------

local robot = {}
local robot_mt = { __index = robot }	-- metatable

local matrix = require 'matrix'

-------------------------------------------------
-- PRIVATE FUNCTIONS
-------------------------------------------------

local function t2r(T)

    local A = matrix {{T[1][1],T[1][2],T[1][3]},
                      {T[2][1],T[2][2],T[2][3]},
                      {T[3][1],T[3][2],T[3][3]}}
    return A

end

local function zeros(m)

    local Z = {}

    for i = 1,m,1 do
        Z[i] = 0
    end

    return Z

end

local function zeros2(m,n)

    local Z = {}

    for i = 1,m,1 do
        Z[i] = {}
        for j = 1,n,1 do
            Z[i][j] = 0
        end
    end

    return Z

end

local function diag(q)

    local n = #q
    local D = {}

    if n == 1 then
        return q
    else
        for i = 1,n,1 do
            D[i] = {}
            for j = 1,n,1 do
                if i == j then
                    D[i][j] = q[i]
                else
                    D[i][j] = 0
                end
            end
        end
    end

    return D

end


-------------------------------------------------
-- PUBLIC FUNCTIONS
-------------------------------------------------

function robot.new(link)	-- constructor

	local newRobot = {
        link = link,
        n = #a
	}

	return setmetatable( newRobot, robot_mt )
end


-------------------------------------------------

function robot:printDH()

    for i = 1,self.n,1 do
	    print(i..': a='..self.link[i].a..
                 ' alpha='..self.link[i].alpha..
                 ' d='..self.link[i].d..
                 ' theta=q'..i)
    end
end

-------------------------------------------------

function robot:gravLoad(q, grav)

    local tau = {}
    local qd = zeros(self.n)
    local qdd = zeros(self.n)

    if grav == nil then
        tau = self:rne(q, qd, qdd)
    else
        tau = self:rne(q, qd, qdd, grav)
    end

    return tau
end

-------------------------------------------------

function robot:coriolis(q, qd)

    local C = {}
    local C2 = matrix {}
    local Csq = {}
    local qdd = zeros(self.n)
    local grav = matrix {{0},{0},{0}}


    -- find the torques that depend on a single finite joint speed,
    -- these are due to the squared (centripetal) terms
    for i = 1,self.n,1 do
        local QD = zeros(self.n)
        QD[i] = 1
        tau = self:rne(q, QD, qdd, grav)
        Csq[i] = tau
        C[i] = matrix {zeros(self.n)}
    end

    -- find the torques that depend on a pair of finite joint speeds,
    -- these are due to the product (Coridolis) terms
    for i = 1,self.n,1 do
        for j = i+1,self.n,1 do
            local QD = zeros(self.n)
            QD[i] = 1
            QD[j] = 1
            tau = matrix {self:rne(q, QD, qdd, grav)}
            temp = 0.5 * (tau - matrix {Csq[j]} - matrix {Csq[i]})
            C[j] = C[j] + temp * qd[i]
            C[i] = C[i] + temp * qd[j]
        end
    end

    C2 = C[1]
    for i = 2,self.n,1 do
       C2 = matrix.concatv(C2, C[i])
    end

    C = matrix.transpose(C2)
    Csq = matrix.transpose(matrix:new(Csq))

    C = C + Csq * matrix:new(diag(qd))

    return C

end
-------------------------------------------------

function robot:rne(q, qd, qdd, grav)

    local vd = {}
    if grav == nil then
        vd = matrix {{0},{0},{9.81}}
    else
        vd = grav
    end

    -- init some variables and compute rotation matrices
    local z0 = matrix {{0},{0},{1}}
    local w  = matrix {{0},{0},{0}}
    local wd = matrix {{0},{0},{0}}

    local link  = {}

    local rm    = {}
    local pm    = {}

    -- init some variables
    for i = 1,self.n,1 do

        link = self.link[i]
        rm[i] = t2r(link:A(q[i]))
        pm[i] = matrix {{link.a},
                       {-link.d*math.sin(link.alpha)},
                       {link.d*math.cos(link.alpha)}}
    end

    local F = {}
    local N = {}

    -- forward interation
    for i = 1,self.n,1 do

        local R  = matrix.transpose(rm[i])
        local P  = pm[i]
        local PC = self.link[i].rc

        local w_  = R*w + z0*qd[i]
        local wd_ = R*wd + matrix.cross(R*w, z0*qd[i]) + z0*qdd[i]
        local vd_ = R*(matrix.cross(wd,P)
                    + matrix.cross(w, matrix.cross(w, P)) + vd)

        -- update variables
        w  = w_
        wd = wd_
        vd = vd_


        local vdc = matrix.cross(wd, PC) +
                    matrix.cross(w, matrix.cross(w, PC)) + vd

        F[i] = vdc*self.link[i].m
        N[i] = self.link[i].I*wd + matrix.cross(w, self.link[i].I*w)

    end

    local f = matrix {{0},{0},{0}}
    local n = matrix {{0},{0},{0}}
    local R = {}
    local P = {}
    local tau = {}

    -- backward iteration
    for i = self.n,1,-1 do

        if i == self.n then
            R = matrix:new(3, 'I')
            P = matrix {{0},{0},{0}}
        else
            R = rm[i+1]
            P = pm[i+1]
        end

        local PC = self.link[i].rc
        local f_ = R*f + F[i]
        local n_ = N[i] + R*n + matrix.cross(PC, F[i]) + matrix.cross(P, R*f)

        f  = f_
        n = n_

        tau[i] = matrix.getelement(n, 3, 1)
    end

    return tau
end

return robot

