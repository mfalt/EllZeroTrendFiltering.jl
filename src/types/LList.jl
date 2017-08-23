type LList{T}
    val::T
    next::LList{T}
    LList{T}() where T = new() #Hack for empty list
    LList{T}(val::T) where T = (l = new(val); l.next = l)
    LList{T}(val::T, l::LList{T}) where T = new(val, l)
end
LList{T}(val::T) = LList{T}(val)
LList{T}(val::T, l::LList{T}) = LList{T}(val, l)

llist{T}(val::T) = LList(val)
llist{T}(val::T, args...) = LList(val, llist(args...))
#Hack for empty set
Base.isempty(l::LList) = !isdefined(l, :next)

function Base.show{T}(io::IO, l::LList{T})
    print(io, "llist(")
    if !isempty(l)
        show(io, l.val)
        if !islast(l)
            for val in l.next
                print(io, ",")
                show(io, val)
            end
        end
    end
    print(io, ")")
end

islast(l::LList) = l === l.next || isempty(l.next)

"""
Insert `val` in `l=llist(v1,v2,...,vn)` so that
`l=llist(v1,val,v2,...,vn)`.
If `l` is empty then `l` becomes `llist(val)`.
"""
function insertafter!{T}(val::T, l::LList{T})
    if isempty(l)
        l.val = val
        l.next = l
    else
        l2 = islast(l) ? LList{T}(val) : LList{T}(val, l.next)
        l.next = l2
    end
    return l
end

Base.start{T}(l::LList{T}) = (l, isempty(l))
Base.done{T}(l::LList{T}, state::Tuple{LList{T},Bool}) = state[2]
Base.next{T}(l::LList{T}, state::Tuple{LList{T},Bool}) = (state[1].val, (state[1].next, islast(state[1])))

function Base.length(l::LList)
    n = 0
    for i in l
        n += 1
    end
    n
end
