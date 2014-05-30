-------------------------------------------------
--
-- robot.lua
--
-------------------------------------------------

local link = {}
local link_mt = { __index = link }	-- metatable

local matrix = require 'matrix'

-------------------------------------------------
-- PRIVATE FUNCTIONS
-------------------------------------------------


-------------------------------------------------
-- PUBLIC FUNCTIONS
-------------------------------------------------

function link.new(a, alpha, d, m, rc, I)	-- constructor

	local newLink = {
        a     = a,
        alpha = alpha,
        d     = d,
        m     = m or 1,
        rc    = rc or matrix {0,0,0},
        I     = I or matrix {{1,0,0},{0,1,0},{0,0,1}}
	}

	return setmetatable( newLink, link_mt )
end


-------------------------------------------------

function link:A(q)

    local sa = math.sin(self.alpha)
    local ca = math.cos(self.alpha);

    local st = math.sin(q)
    local ct = math.cos(q);

    local a  = self.a;
    local d  = self.d;

    local T = matrix {{ct, -st, 0, a},
                      {st*ca, ct*ca, -sa, -sa*d},
                      {st*sa, ct*sa, ca, ca*d},
                      {0, 0, 0, 1}}

    return T
end

-------------------------------------------------

function link:printDH()

	    print('a='..self.a..
              ' alpha='..self.alpha..
              ' d='..self.d..
              ' theta=q')
end

return link
