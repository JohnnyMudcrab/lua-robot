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

function robot:rne(q, Qd, Qdd)

    -- set velocity and acceleration to zero if not given
    local qd = {}
    local qdd = {}

    if Qdd == nil then
        for i = 1,self.n,1 do qdd[i] = 0 end
    else
        qdd = Qdd
    end

    if Qd == nil then
        for i = 1,self.n,1 do qd[i] = 0 end
    else
        qd = Qd
    end

    -- init some variables and compute rotation matrices
    local z0 = matrix {{1},{2},{3}}
    local w  = matrix {{0},{0},{0}}
    local wd = matrix {{0},{0},{0}}
    local vd = matrix {{0},{0},{9.81}}

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
        local wd_ = R*wd
        local vd_ = R*(matrix.cross(wd,P)
                    + matrix.cross(w, matrix.cross(w, P)) + vd)

        -- update variables
        w  = w_
        wd = wd_
        vd = vd_

        local vdc = matrix.cross(wd, PC) +
                    matrix.cross(w, matrix.cross(w, PC)) + vd

        F[i] = vdc*self.link[i].m
        N[i] = self.link[i].I*wd + matrix.cross(w, self.link[i].I)

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
        print(tau[i])
    end

    return tau
end

return robot

