program mc_nvt
    implicit none
    integer, parameter :: numPart = 1000
    real, parameter :: density = 0.844
    real :: disp
    real :: init_energy
    integer :: steps


    real, dimension(1:3, 1:numPart) :: coords

    call init_coords(coords)
    !call print_coords(coords)
    call run_simulation(coords)

    contains
        subroutine init_coords(coords)
            !Hard coded, next revision make this more robust
            implicit none
            !specifically for numPart = 100
            !integer, parameter :: nCube = 5

            integer, parameter :: nCube = 10

            integer :: i
            real    :: L
            real, dimension(1:3) :: indexer = (/ 0.0, 0.0, 0.0 /)
            real, dimension(1:3) :: offset = (/ 0.5, 0.5, 0.5 /)
            real, dimension(1:3, 1:numPart) :: coords
            real :: init_energy

            L = (numPart / density)**(1.0 / 3.0)

            do i = 1, numPart
                coords(:, i) = (indexer + offset) * (L / nCube)
                !advance the indexer
                indexer(1) = indexer(1) + 1
                if (indexer(1) == nCube) then
                    indexer(1) = 0.0;
                    indexer(2) = indexer(2) + 1
                    if(indexer(2) == nCube) then
                        indexer(2) = 0.0
                        indexer(3) = indexer(3) + 1
                    end if
                end if
            end do

        end subroutine

        !Prints out the current xyz coordinates of all atoms
        subroutine print_coords(coords)
            implicit none
            real, dimension(1:3, 1:numPart) :: coords
            integer :: i,j

            do j = 1, numPart
                write (*, "(a)", advance = "no") "Coordinate : "
                write (*, "(I4)") j
                do i = 1, 2
                    write (*, "(f10.5)", advance = "no") coords(i, j)
                end do
                write (*, "(f10.5)", advance = "yes") coords(3, j)
            end do

        end subroutine

        subroutine print_xyz(coords)
            implicit none
            real, dimension(1:3, 1:numPart) :: coords
            integer :: i, j

            write (*, *) numPart
            do j = 1, numPart
                write (*, "(a)", advance = "no") "C"
                do i = 1, 2
                    write (*, "(f10.5)", advance = "no") coords(i, j)
                end do
                write (*, "(f10.5)", advance = "yes") coords(3, j)
            end do

        end subroutine

        !This takes c and adjusts it to be within 0 to L for all axes
        subroutine adjust_PBC(c, coords)
            implicit none
            real, dimension(1:3) :: c
            real, dimension(1:3, 1:numPart) :: coords
            integer :: i
            real :: L

            L = (numPart / density)**(1.0 / 3.0)

            do i = 1, 3
                do while (c(i) < 0)
                    c(i) = c(i) + L
                end do
                c(i) = mod(c(i), L)
            end do

        end subroutine

        !returns distance as a vector with x, y, z coordinates
        !If the distance is over L, this means there is another version of the same
        !Particle(s) that is closer, and we put take that smaller distance
        function distance(c)
            implicit none
            real, dimension(1:3) :: c
            real, dimension(1:3) :: distance
            integer :: i
            real L, hL

            L = (numPart / density)**(1.0/3.0)
            hL = L / 2

            distance = c
            do i = 1, 3
                if (distance(i) > hL) then
                    distance(i) = distance(i) - L
                else if (distance(i) < -hL) then
                    distance(i) = distance(i) + L
                end if
            end do
            return
        end function

        !L-J 6 12
        real function energy(coords)
            implicit none
            real, dimension(1:3, 1:numPart) :: coords
            real :: L, dist_sq
            real, dimension(1:3) :: dist
            integer :: i,j

            energy = 0.0
            L = (numPart / density)**(1.0 / 3.0)
            !all pairs
            do i = 1, (numPart - 1)
                do j = (i + 1), numPart
                    dist = coords(:, i) - coords(:, j)
                    dist = distance(dist)
                    dist_sq = dist(1)**2 + dist(2)**2 + dist(3)**2

                    energy = energy + ((1 / dist_sq)**6 - (1 / dist_sq)**3)
                end do
            end do

            energy = energy * 4
        return
        end function

        !loc is index of the coordinate to be changed
        real function energy_delta(coords, loc, new_coord)
            implicit none
            real, dimension(1:3, 1:numPart) :: coords
            integer :: loc
            real, dimension(1:3) :: new_coord

            real, dimension(1:3) :: old_coord
            integer :: i
            real :: L, new_dist_sq, old_dist_sq, energy_old, energy_new
            real, dimension(1:3) :: new_dist, old_dist


            old_coord = coords(:, loc)

            energy_old = 0.0
            energy_new = 0.0
            L = (numPart / density) ** (1.0 / 3.0)

            !only pairs with the particle of interest
            do i = 1, numPart
                if (i /= loc) then
                    new_dist = coords(:, i) - new_coord
                    old_dist = coords(:, i) - old_coord

                    new_dist = distance(new_dist)
                    old_dist = distance(old_dist)

                    new_dist_sq = new_dist(1)**2 + new_dist(2)**2 + new_dist(3)**2
                    old_dist_sq = old_dist(1)**2 + old_dist(2)**2 + old_dist(3)**2

                    energy_new = energy_new + ((1 / new_dist_sq)**6 - (1 / new_dist_sq)**3)
                    energy_old = energy_old + ((1 / old_dist_sq)**6 - (1 / old_dist_sq)**3)
                end if
            end do

            energy_delta = energy_new - energy_old
            energy_delta = energy_delta * 4
            return
        end function

        subroutine radial_dist_sample(g_r, nhis, del_g, L, coords)
            implicit none
            real, dimension(0:(nhis - 1)) :: g_r

            integer :: nhis
            real :: del_g, L
            real, dimension(1:3, 1:numPart) :: coords
            real, dimension(1:3) :: dist_v
            real :: dist_s
            integer :: bin_index, i, j

            do i = 1, (numPart - 1)
                do j = (i + 1), numPart
                    dist_v = coords(:, i) - coords(:, j)
                    dist_v = distance(dist_v) !adjust for Periodic Boundary conditions
                    dist_s = sqrt(dist_v(1)**2 + dist_v(2)**2 + dist_v(3)**2)
                    if (dist_s <= (L / 2)) then
                        bin_index = int(dist_s / del_g)
                        g_r(bin_index) = g_r(bin_index) + 2 !count both in the pair
                    end if
                end do
            end do

        end subroutine

        !modifies g_r to normalize it
        !nPart =
        subroutine radial_dist_normalize(g_r, nhis, del_g, numPart, ngr)
            implicit none
            real, dimension(0:(nhis - 1)) :: g_r
            integer :: nhis, numPart, ngr
            real :: del_g

            real :: vb, nid!, r
            integer :: i
            real, parameter :: pi = 3.141592

            !determine g(r)
            do i = 0, (nhis - 1)
                !r = del_g * (i + 0.5) !average distance of a histogram bin
                vb = (((i+1)**3) - (i**3)) * (del_g**3) !volume of this bin
                !nid = (4/3) * pi * vb * density !???"::
                nid = 4 * pi * (del_g * (i + .5)**2) * (del_g)
                !nid = pi * vb * density
                g_r(i) = g_r(i) / (ngr * numPart * nid) !normalize g_r
            end do
        end subroutine

        subroutine print_radial_dist(g_r, nhis, del_g)
            implicit none
            real, dimension(0:(nhis - 1)) :: g_r
            integer :: nhis
            real :: del_g

            integer :: i
            real :: r

            write (*, "(a)", advance = "no") "          r          "
            write (*, *) "g(r)"
            do i = 0, (nhis - 1)
                r = del_g * (i + 0.5)
                write (*, "(f15.8)", advance = "no") r
                write (*, "(a)", advance = "no") " "
                write (*, *) g_r(i)
            end do

        end subroutine

        subroutine run_simulation(coords)
            implicit none
            real, dimension(1:3, 1:numPart) :: coords
            real, parameter :: temp = 0.71
            real, parameter :: beta = 1.0 / temp
            real :: disp
            real :: energy_system
            real :: delta_energy
            integer, parameter :: steps = 4000
            integer :: i,j

            real, parameter :: maxDisp = 0.1
            real, dimension(1:3) :: new_direction
            real, dimension(1:3) :: new_coord
            integer :: counter = 0
            integer, parameter :: seed = 36459!32415

            !g(r) parameters

            integer :: ngr = 0
            integer, parameter ::nhis = 1000
            real :: del_g, L
            real, dimension(0:(nhis - 1)) :: g_r

            L = (numPart / density) ** (1.0 / 3.0)

            del_g = (L * sqrt(3.0)) /  nhis
            do i = 0, (nhis - 1)
                g_r(i) = 0
            end do


            energy_system = energy(coords)
            !write (*, *) "Initial Energy:"
            !write (*, *) energy_system

            call srand(seed)
            do i = 1, steps
                !if ((mod(i, 3998)) == 0) then
                    !write (*, *) i
                !    call radial_dist_sample(g_r, nhis, del_g, L, coords)
                !    ngr = ngr + 1
                !end if
                !if ((mod(i, 1000)) == 0) then
                !    write (*, *) i
                !end if

                do j = 1, numPart
                    new_direction = (/ rand() - .5, rand() - .5, rand() - .5 /)
                    new_coord = coords(:, j) + maxDisp * new_direction

                    call adjust_PBC(new_coord, coords)
                    delta_energy = energy_delta(coords, j, new_coord)

                    if (rand() < exp(-1* beta * delta_energy)) then
                        coords(:, j) = new_coord
                        energy_system = energy_system + delta_energy
                        counter = counter + 1
                    end if
                end do
            end do
            !write (*, *) "Final Energy:"
            !write (*, *) energy_system
            call radial_dist_normalize(g_r, nhis, del_g, numPart, ngr)
            call print_radial_dist(g_r, nhis, del_g)
            !write (*, *) "Done"
            call print_xyz(coords)
        end subroutine

end program mc_nvt
