using Delaunay
using Plots
using Statistics
using IterTools

function lenOfVec(x1::Float64,x2::Float64)
    return sqrt(x1*x1 + x2*x2)
end

function plotMesh(tri)
    forPlotx = []
    forPloty = []

    for i in eachrow(tri.simplices)
        push!(forPlotx, points[i[1],1])
        push!(forPloty, points[i[1],2])
        push!(forPlotx, points[i[2],1])
        push!(forPloty, points[i[2],2])
        push!(forPlotx, points[i[3],1])
        push!(forPloty, points[i[3],2])
        push!(forPlotx, points[i[1],1])
        push!(forPloty, points[i[1],2])
        push!(forPlotx, missing)
        push!(forPloty, missing)
    end

    plot!(forPlotx,forPloty)
end

function generatePoints(x_start::Float64,x_stop::Float64,y_start::Float64,y_stop::Float64,step::Float64)
    x_points = [x for x in x_start:step:x_stop]
    y_points = [y for y in y_start:step:y_stop]

    points = zeros(length(x_points)*length(y_points),2)

    indx = 1
    for x in x_points
        for y in y_points
            points[indx, 1] = x
            points[indx, 2] = y
            indx += 1
        end
    end

    return points
end

function isInConvexHull(edge1::Int64,edge2::Int64, tri::Triangulation)
    for (i,j) in eachrow(tri.convex_hull)
        if i == edge1 && j == edge2
            return true
        end
    end
    
    return false
end

function createVecFrom2Points(a, b)
    return [b[1]-a[1],b[2]-a[2]]
end

function vecDot(a,b)
    return a[1]*b[1] + a[2]*b[2]
end

mutable struct EdgeC
    edgeS1::Int64
    edgeS2::Int64
    edgeP1::Array{Float64}
    edgeP2::Array{Float64}
    borderEdge
    gDiff_f::Float64
    Tfvec
    neighbor::Int64
    S_fvec
    g_C::Float64
    g_f::Float64
    EdgeC(edgeS1,edgeS2,edgeP1,edgeP2) = new(edgeS1,edgeS2,edgeP1,edgeP2)
end

mutable struct NodeC
    p::Array{Float64}
    V_C::Float64
    edges::Array{EdgeC}
    NodeC(p,V_C,edges) = new(p, V_C,edges)
end

function computeGrad(V_C, indFlux, flux, node)
    gradient = [0.,0.]

    for edge in node.edges
        if(edge.borderEdge)
            phi_f = boundaryCondition
            gradient[1] += phi_f * edge.S_fvec[1]
            gradient[2] += phi_f * edge.S_fvec[2]
        else
            phi_f = edge.g_C * flux[indFlux] + (1-edge.g_C)*flux[edge.neighbor]
            gradient[1] += phi_f * edge.S_fvec[1]
            gradient[2] += phi_f * edge.S_fvec[2]
        end
    end
    #println([gradient[1] / V_C, gradient[2] / V_C])
    return [gradient[1] / V_C, gradient[2] / V_C]
    #gradPhi_C = 1 / V_C * sum phi_f * S_fvec
end

function generateNodes(tri::Triangulation, points::Array{Float64,2})::Array{NodeC,1}
    nodesC = NodeC[]
    for (ind,s) in enumerate(eachrow(tri.simplices))
        p1 = mean(points[s,1])
        p2 = mean(points[s,2])
        V_C = 1/2 * 
                (points[s[1],1] * (points[s[2],2]-points[s[3],2]) +
                points[s[2],1] * (points[s[3],2]-points[s[1],2]) + 
                points[s[3],1] * (points[s[1],2]-points[s[2],2]))
        newNode = NodeC([p1,p2],V_C,EdgeC[])
        push!(nodesC, newNode)
    end
    
    for (ind,s) in enumerate(eachrow(tri.simplices))

        for (edgeS1,edgeS2) in subsets(s,2)
            edgeP1 = points[edgeS1,:]
            edgeP2 = points[edgeS2,:]
            newEdge = EdgeC(edgeS1,edgeS2,edgeP1,edgeP2)
            
            S_fvec = [edgeP1[2]-edgeP2[2], -(edgeP1[1]-edgeP2[1])]
            edgeCenter = [(edgeP1[1]+edgeP2[1]) / 2,(edgeP1[2]+edgeP2[2]) / 2]

            d_Cfvec = createVecFrom2Points(nodesC[ind].p[:],edgeCenter[:])
            d_Cf = lenOfVec(d_Cfvec[1], d_Cfvec[2])
            
            edgeCenter_dot_d_Cf = vecDot(edgeCenter,d_Cfvec)
            if edgeCenter_dot_d_Cf < 0
                S_fvec[1] = -S_fvec[1]
                S_fvec[2] = -S_fvec[2]
            end
            S_f = lenOfVec(S_fvec[1],S_fvec[2])
            nvec = [S_fvec[1] / S_f, S_fvec[2] / S_f]

            borderEdge = false
            #edgeS1, edgeS2 in tri.convex_hull -> edgeS1,edgeS2 jsou postranni
            #pokud tam nejsou, tak zkusit neighbors a kde jsou ty stejny
            if (isInConvexHull(edgeS1,edgeS2,tri))
                borderEdge = true

                evec = [d_Cfvec[1] / d_Cf, d_Cfvec[2] / d_Cf]
                Efvec = [S_f * evec[1], S_f * evec[2]]
                Ef = lenOfVec(Efvec[1],Efvec[2])

                Tfvec = [(nvec[1] - evec[1])*S_f, (nvec[2] - evec[2])*S_f]
                Tf = lenOfVec(Tfvec[1], Tfvec[2])

                gDiff_f = Ef / d_Cf
                newEdge.gDiff_f = gDiff_f
                newEdge.Tfvec = Tfvec
            else
                for neighbor in tri.neighbors[ind,:]
                    if (neighbor!=0)
                        simplOfNei = tri.simplices[neighbor,:]
                        if edgeS1 in simplOfNei && edgeS2 in simplOfNei
                            d_CFvec = createVecFrom2Points(nodesC[ind].p[:],nodesC[neighbor].p[:])
                            d_CF = lenOfVec(d_CFvec[1],d_CFvec[2])
                            d_fFvec = createVecFrom2Points(edgeCenter[:], nodesC[neighbor].p[:])
                            d_fF = lenOfVec(d_fFvec[1], d_fFvec[2])
                            
                            evec = [d_CFvec[1] / d_CF, d_CFvec[2] / d_CF]
                            Efvec = [S_f * evec[1], S_f * evec[2]]
                            Ef = lenOfVec(Efvec[1],Efvec[2])
                            
                            Tfvec = [(nvec[1] - evec[1])*S_f, (nvec[2] - evec[2])*S_f]
                            Tf = lenOfVec(Tfvec[1], Tfvec[2])

                            g_f = d_Cf / (d_Cf + d_fF)
                            g_C = d_fF / d_CF
                            gDiff_f = Ef/d_CF

                            newEdge.gDiff_f = gDiff_f
                            newEdge.Tfvec = Tfvec
                            newEdge.neighbor = neighbor
                            newEdge.g_C = g_C
                            newEdge.g_f = g_f
                        end
                    end
                end
            end
            newEdge.borderEdge = borderEdge
            newEdge.S_fvec = S_fvec
            #println(newEdge)
            #println(nodesC[ind].edges)
            push!(nodesC[ind].edges, newEdge)
        end
        
    end
    
    return nodesC
end

function computeCoef!(flux, nodes)

    for (ind_node,node) in enumerate(nodes)

        #fluxC_copy = flux[ind_node]

        a_f = 0.0
        a_C = 0.0
        b_C = Sc * node.V_C
    
        gradC = computeGrad(node.V_C, ind_node, flux, node)
        #println(gradC, "gradC")
        for ed in node.edges
            if ed.borderEdge
                #gDiff_b = Eb / d_Cb
                #=
                flux[ind_node] = boundaryCondition
                break
                =#
                FluxC_b = Gamma * ed.gDiff_f
                FluxV_b = -Gamma * ed.gDiff_f * boundaryCondition - Gamma * vecDot(computeGrad(node.V_C, ind_node, flux, node), ed.Tfvec)
                #println(FluxC_b, " FluxC_b ", FluxV_b, " FluxV_b")
                #a_f = FluxF_f
                a_C += FluxC_b
                b_C += -FluxV_b

            else
                #println(ind_node," ", ed.neighbor)
                a_f += -Gamma * ed.gDiff_f * flux[ed.neighbor]
                a_C += Gamma * ed.gDiff_f
    
                gradF = computeGrad(nodes[ed.neighbor].V_C, ed.neighbor, flux, nodes[ed.neighbor])
                gradPhi_f = [ed.g_f * gradF[1] + (1-ed.g_f)*gradC[1], ed.g_f * gradF[2] + (1-ed.g_f)*gradC[2]]
                
                b_C += Gamma * vecDot(gradPhi_f, ed.Tfvec)
                #println(gradF, " gradF ")
            end
        end
        
        #= TODO: Some Underrelaxation - should be this easy
        a_C = a_C / lambda
        b_C = b_C + ((1-lambda)*a_C)*fluxC_copy
        =#
        #println(a_C, " a_C ", a_f, " a_f ", b_C, "b_C")
        flux[ind_node] = (-a_f + b_C) / (a_C) 
    
    end
end

# ------- PARAMETERS -----------
const Sc = 0.          # Linearization of source term
const Sp = 0.
const Gamma = 0.02      # Assuming equal Gamma (diffusive coefficient) throughout the domain
const rho = 1.          # Assuming equal density throughout the domain
const vel_u = 1.        # Velocity of flow along X direction
const vel_v = 4.        # Velocity of flow along Y direction
const boundaryCondition = 2.   # boundaryCondition - for now const
#const lambda = 0.2     # Underrelaxation const

points = generatePoints(-1.,1.,-1.,1.,0.1)
tri = delaunay(points)

nodes = generateNodes(tri, points)

flux = fill(1., size(nodes)[1])
for i in 1:5
    computeCoef!(flux,nodes)
end

#plotMesh(tri)
#scatter!([i.p[1] for i in nodes], [i.p[2] for i in nodes])
#contour([i.p[1] for i in nodes], [i.p[2] for i in nodes], [fl for fl in flux])
#plot!(points[:,1], points[:,2])
#plot!([i for i in 1:size(nodes)[1]], [fl for fl in flux])
#l = @layout [a{0.7w} b; c{0.2h}]
#plot([i.p[1] for i in nodes], [i.p[2] for i in nodes], [fl for fl in flux], st = [:surface, :contourf], layout = l)
surface([i.p[1] for i in nodes],[i.p[2] for i in nodes],[fl for fl in flux])