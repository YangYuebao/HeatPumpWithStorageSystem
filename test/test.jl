struct mystruct
    a::Int
    b::Int
end
ms = mystruct(1,2)
macro generate_a_b(ms)
    ex = quote
        a=$(ms).a
        b=$(ms).b
    end
    return esc(ex) 
end
a=0
b=0
@generate_a_b(ms)
