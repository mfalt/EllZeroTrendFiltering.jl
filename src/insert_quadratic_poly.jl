
"""
    insert_quadratic!(Λ::PiecewiseQuadratic{T}, μ::QuadraticPolynomial{T}) where {T}
Inserts a quadratic polynomial μ into the linked list `Λ`, which represents a piecewise
quadratic polynomial, so that the new `Λ` satisfies
`Λ(y) := min{ Λ(y) ,μ(y) } ∀ y`
"""
function insert_quadratic!(Λ::PiecewiseQuadratic{T}, μ::QuadraticPolynomial{T}) where T

    if Λ.next.left_endpoint == Inf # I.e. the piecewise quadratic object is empty, perhaps better to add dummy polynomial
        μ.has_been_used = true
        insert(Λ, μ, -Inf) # FIXME: Is specific insert function needed?
        return
    end

    λ_prev = Λ # Points to the list head which is just a node with NaN data
    λ_curr = Λ.next # The leftmost segment of the list

    while λ_curr.left_endpoint != Inf #???
        #global counter2 += 1
        DEBUG && println(Λ)

        left_endpoint = λ_curr.left_endpoint
        right_endpoint = get_right_endpoint(λ_curr)


        Δa = μ.a - λ_curr.p.a
        Δb = μ.b - λ_curr.p.b
        Δc = μ.c - λ_curr.p.c

        b2_minus_4ac =  Δb^2 - 4*Δa*Δc

        if Δa > 0 # μ has greater curvature, i.e., μ is smallest in the middle if intersect
            if b2_minus_4ac <= ACCURACY
                # Zero (or one) intersections, old quadratic is smallest, just step forward
                DEBUG && println("No intersections, old quadratic is smallest, Δa > 0, breaking.")
                break
            else

                # Compute the intersections
                term1 = -(Δb / 2 / Δa)
                term2 = sqrt(b2_minus_4ac) / 2 / Δa # Δa > 0
                root1, root2 = term1-term2, term1+term2

                DEBUG && println("Δa > 0   root1:", root1, "   root2:", root2)

                # Check where the intersections are and act accordingly
                if root1 >= right_endpoint
                    DEBUG && println("Two intersections to the right: old quadratic is smallest")
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                elseif root2 <= left_endpoint
                    DEBUG && println("Two intersections to the left: old quadratic is smallest")
                    break # There will be no more intersections since Δa > 0
                elseif root1 <= left_endpoint && root2 >= right_endpoint
                    DEBUG && println("One intersections on either side, new quadratic is smallest")
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                elseif root1 > left_endpoint && root2 < right_endpoint
                    DEBUG && println("Two intersections: old-new-old")
                    λ_prev, λ_curr = update_segment_old_new_old(λ_curr, μ, root1, root2)
                    break # There will be no more intersections since Δa > 0
                elseif root1 > left_endpoint
                    DEBUG && println("Root 1 within the interval: old-new")
                    λ_prev, λ_curr = update_segment_old_new(λ_curr, μ, root1)
                elseif root2 < right_endpoint
                    DEBUG && println("Root 2 within the interval: new-old")
                    λ_prev, λ_curr = update_segment_new_old(λ_prev, λ_curr, μ, root2)
                    break # There will be no more intersections since Δa > 0
                else
                    error("Shouldn't end up here")
                end
            end

        elseif Δa < 0 # μ has lower curvature, i.e., μ is smallest on the sides
            if b2_minus_4ac <= ACCURACY
                # Zero (or one) roots
                DEBUG && println("Δa < 0   b2_minus_4ac: ", b2_minus_4ac)
                λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
            else
                # Compute the intersections
                term1 = -(Δb / 2 / Δa)
                term2 = sqrt(b2_minus_4ac) / 2 / Δa # < 0
                root1, root2 = term1+term2, term1-term2
                DEBUG && println("Δa < 0   root1:", root1, "   root2:", root2)

                # Check where the intersections are and act accordingly
                if root1 >= right_endpoint || root2 <= left_endpoint
                    # No intersections, μ is smallest
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                elseif root1 <= left_endpoint && root2 >= right_endpoint
                    # No intersections, old quadratic is smallest, just step forward
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                elseif root1 > left_endpoint && root2 < right_endpoint
                    # Two intersections within the interval
                    λ_prev, λ_curr = update_segment_new_old_new(λ_prev, λ_curr, μ, root1, root2)
                elseif root1 > left_endpoint
                    λ_prev, λ_curr = update_segment_new_old(λ_prev, λ_curr, μ, root1)
                elseif root2 < right_endpoint
                    λ_prev, λ_curr = update_segment_old_new(λ_curr, μ, root2)
                else
                    error("Shouldn't end up here")
                end
            end
        else # Δa == 0.0
            DEBUG && println("Δa == 0")
            DEBUG2 && println("μ       : $μ")
            DEBUG2 && println("Δb      : $Δb")
            DEBUG2 && println("Δc      : $Δc")
            if Δb == 0
                if Δc >= 0
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                else
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                end
                continue
            end

            root = -Δc / Δb
            DEBUG && println("root    : $root")
            DEBUG && println("left_endpoint : $left_endpoint")
            DEBUG && println("right_endpoint : $right_endpoint")
            # TODO think through this more and reduce to one case?
            # These two cases only seem to happen for LTI systems
            if root < -1e10 # Don't allow intersections far away for numerical reasons
                DEBUG && println("root    : $root < -1e10, Special case 1")
                if Δc >= 0 # intersection to left we and new is worse shouldnt have to do anything
                    # TODO Break here?
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                else # New is better, keep new
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                end
            elseif root > 1e10
                DEBUG && println("root    : $root > 1e10, Special case 2")
                if Δc >= 0 # intersection to RIGHT we and new is better far to the right, ignore this?
                    # TODO Something else?
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                else # New is better, keep new
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                end
            else
                if Δb > 0
                    if root < left_endpoint
                        λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                    elseif root > right_endpoint
                        λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                    else
                        λ_prev, λ_curr = update_segment_new_old(λ_prev, λ_curr, μ, root)
                    end
                else
                    if root < left_endpoint
                        λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                    elseif root > right_endpoint
                        λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                    else
                        DEBUG2 && println("Special case")
                        λ_prev, λ_curr = update_segment_old_new(λ_curr, μ, root)
                    end
                end
            end
        end

    end
    return
end

@inline function update_segment_do_nothing(λ_curr)
    return λ_curr, λ_curr.next
end

@inline function update_segment_new(λ_prev, λ_curr, μ)
    μ.has_been_used = true # TODO: Could be moved to the second clause
    if λ_prev.p === μ
        λ_prev.next = λ_curr.next
        v1, v2 = λ_prev, λ_curr.next
    else
        λ_curr.p = μ
        v1, v2 = λ_curr, λ_curr.next
    end
    return v1, v2
end

@inline function update_segment_old_new(λ_curr, μ, root)
    # ... λ_curr | λ_new ...
    μ.has_been_used = true
    λ_new =  PiecewiseQuadratic(μ, root, λ_curr.next)
    λ_curr.next = λ_new
    return λ_new, λ_new.next
end

@inline function update_segment_new_old(λ_prev, λ_curr, μ, root)
    # ... λ_prev | (new segment) | λ_curr ...
    if λ_prev.p === μ
        λ_curr.left_endpoint = root
    else
        μ.has_been_used = true
        λ_prev.next = PiecewiseQuadratic(μ, λ_curr.left_endpoint, λ_curr)
        λ_curr.left_endpoint = root
    end
    return λ_curr, λ_curr.next
end

@inline function update_segment_new_old_new(λ_prev, λ_curr, μ, root1, root2)
    # Insert new segments with μ before and after λ_curr
    # ... λ_prev : (new segment) | λ_curr | (new segment) ...
    # ( μ.has_been_used = true is set inside the called funcitons )
    update_segment_new_old(λ_prev, λ_curr, μ, root1)
    return update_segment_old_new(λ_curr, μ, root2)
end

@inline function update_segment_old_new_old(λ_curr, μ, root1, root2)
    # ... λ_curr | λ_1  | λ_2 ...
    # The second new piece contains a copy of the polynomial in λ_curr

    μ.has_been_used = true
    λ_2 =  PiecewiseQuadratic(λ_curr.p, root2, λ_curr.next)
    λ_1 =  PiecewiseQuadratic(μ, root1, λ_2)
    λ_curr.next = λ_1
    return λ_2, λ_2.next
end


"""
Given a piecewise quadratic V_Λ represented by Λ, a polynomial μ, and a real number ζ
this function evaluates if
μ(y) > V_Λ(y) + ζ  ∀ y
""" # FIXME: Group polynomial and constant in tuple?
function poly_minus_constant_is_greater(Λ::PiecewiseQuadratic{T}, μ::QuadraticPolynomial{T}, ζ::Real) where T

    if Λ.next.left_endpoint == Inf
        return false
    end

    for λ_curr in Λ

        left_endpoint = λ_curr.left_endpoint
        right_endpoint = get_right_endpoint(λ_curr)

        Δa = μ.a - λ_curr.p.a
        Δb = μ.b - λ_curr.p.b
        Δc = (μ.c - ζ) - λ_curr.p.c

        b2_minus_4ac =  Δb^2 - 4*Δa*Δc

        if Δa > 0 # μ has greater curvature, i.e., μ is smallest in the middle if intersect
            #println("Δa > 0")
            if b2_minus_4ac <= 0
                # Zero (or one) intersections, old quadratic is smallest, just step forward
                #println("No intersections, old quadratic is smallest, Δa > 0, breaking.")
                return true
            else

                # Compute the intersections
                term1 = -(Δb / 2 / Δa)
                term2 = sqrt(b2_minus_4ac) / 2 / Δa # Δa > 0
                root1, root2 = term1-term2, term1+term2

                # Check where the intersections are and act accordingly
                if root1 >= right_endpoint
                    continue
                elseif root2 <= left_endpoint
                    return true # There will be no more intersections since Δa > 0
                else
                    return false
                end
            end

        elseif Δa < 0 # μ has lower curvature, i.e., μ is smallest on the sides
            #println("Δa < 0")
            if b2_minus_4ac <= 0
                # Zero (or one) roots
                return false
            else
                # Compute the intersections
                term1 = -(Δb / 2 / Δa)
                term2 = sqrt(b2_minus_4ac) / 2 / Δa # < 0
                root1, root2 = term1+term2, term1-term2

                # Check where the intersections are and act accordingly
                if root1 <= left_endpoint && root2 >= right_endpoint
                    # No intersection on either side of the interval, and μ is smallest
                    continue
                else
                    return false
                end
            end
        else # Δa == 0.0
            DEBUG2 && println("Δa == 0 : $μ")
            if Δb == 0
                if Δc >= 0

                else
                    return false
                end
                continue
            end

            root = -Δc / Δb
            if Δb > 0
                if root < left_endpoint

                elseif root > right_endpoint
                    return false
                else
                    return false
                end
            else
                if root < left_endpoint
                    return false
                elseif root > right_endpoint

                else
                    DEBUG2 && println("Special case")
                    return false
                end
            end
        end

    end
    return true
end
