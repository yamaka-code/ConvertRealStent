program main
    implicit none


    integer, Parameter :: nodeNum = 50 !ここをプロット点数に合わせて変更

    integer :: i, j, k, ltemp, count, number
    integer, Parameter :: addPNum = 3 !プロットした点を補間する点の数
    integer, Parameter :: addWNum = 3
    integer, Parameter :: inputWNum = 6 !インプットするwireの本数
    integer, Parameter :: totalWNum = inputWNum * (addWNum + 1) 
    integer, Parameter :: allPNum =  2 + (nodeNum - 2)*(addPNum + 1) - addPNum 
 
    double precision, allocatable, dimension(:,:,:) :: pointLeft, pointRight
    double precision, dimension(1:3) :: center = 0.0d0, adjust = 0.0d0
    double precision :: distance = 0.0d0

    character filename*128, test*100

    allocate(pointLeft(1:3,1:nodeNum,1:totalWNum),pointRight(1:3,1:nodeNum,1:totalWNum))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  input
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, inputWNum
        write(filename, '("LandmarksLeft", I1.1, ".landmarkAscii")') i
        open(10, file=filename, status='old')

        do j = 1, 14
            read (10,'()')        
        enddo

        do j = 1, nodeNum
            read(10, *) pointLeft(1:3,j,(i-1)*(addWNum+1)+1) 
            write(*,*) pointLeft(1:3,j,(i-1)*(addWNum+1)+1) 
        enddo
        close(10)
    enddo    

    do i = 1, inputWNum
        write(filename, '("LandmarksLeft", I1.1, ".landmarkAscii")') i
        open(11, file=filename, status='old')

        do j = 1, 14
            read (11,'()')        
        enddo

        do j = 1, nodeNum
            read(11, *) pointRight(1:3,j,(i-1)*(addWNum+1)+1) 
        !    write(*,*) pointRight(1:3,j,(i-1)*(addWNum+1)+1) 
        enddo
        close(11)
    enddo  

    !!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  中心点とその位置での半径（最小値、最大値、平均値）を計算（交点のみ）
    !
    !!!!!!!!!!!!!!!!!!!!!!!!

    open(13,file = 'RadiusDistributionIntersection.csv')
    write(13,*)'center(1),center(2),center(3),Distance,RadiiMin,RadiiMax,RadiiAve'

    do i = 1, nodeNum
        center = 0.0d0
        do j = 1,inputWNum
            center(:) = center(:) + pointLeft(:, i, (j-1)*(addWNum + 1) +1)
        enddo

        do j = 1, inputWNum
            center(:) = center(:) + pointRight(:, i, (j-1)*(addWNum + 1) +1)
        enddo

        center(:) = center(:) / dble(inputWNum * 2)

        radiiMin = 1000000.0d0
        radiiMax = 0.0d0
        radiiAve = 0.0d0

        do j = 1, inputWNum
            temp = sqrt((pointLeft(1,i,(j-1)*(addWNum+1)+1)-center(1))**2+ &
                &(pointLeft(2,i,(j-1)*(addWNum+1)+1)-center(2))**2+ &
                &(pointLeft(3,i,(j-1)*(addWNum+1)+1)-center(3))**2)
            
            if(radiiMin > temp) radiiMin = temp
            if(radiiMax < temp) radiiMax = temp
            radiiAve = radiiAve + temp
        enddo

        do j = 1, inputWNum
            temp = sqrt((pointRight(1,i,(j-1)*(addWNum+1)+1)-center(1))**2+ &
                &(pointRight(2,i,(j-1)*(addWNum+1)+1)-center(2))**2+ &
                &(pointRight(3,i,(j-1)*(addWNum+1)+1)-center(3))**2)
            
            if(radiiMin > temp)radiiMin = temp
            if(radiiMax < temp)radiiMax = temp
            radiiAve = radiiAve + temp
        enddo

        radiiAve = radiiAve / dble(inputWNum * 2) 

        if(i.ne.1)then
            distance = distance + sqrt((center(1)-adjst(1))**2 + (center(2)-adjst(2))**2 + (center(3)-adjst(3))**2)
        endif
        adjst(:) = center(:)
        write(13,*) center(1),',',center(2),',',center(3),',',distance,',',radiiMin,',',radiiMax,',',radiiAve

    enddo
    close(13)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  tip以外でRightのLeftに近い点をLeftに合わせる
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 3,nodeNum-2
        do j = 1, inputWNum
            distance = 10000000.0d0
            do k = 1, inputWNum
                if(pointRight(1, i, (j - 1)*(interpNum + 1) + 1) == 0.0d0) write(*,*)"zero!"
                if(pointLeft(1, i, (k - 1)*(interpNum + 1) + 1) == 0.0d0) write(*,*)"zero!"
                
                temp = sqrt((pointRight(1,i,(j-1)*(interpNum+1)+1) - pointLeft(1,i,(k-1)*(interpNum+1)+1))**2 + &
                    &(pointRight(2,i,(j-1)*(interpNum+1)+1) - pointLeft(2,i,(k-1)*(interpNum+1)+1))**2 + &
                    &(pointRight(3,i,(j-1)*(interpNum+1)+1) - pointLeft(3,i,(k-1)*(interpNum+1)+1))**2 )
                
                if(temp < distance)then
                    distance = temp
                    number = k
                endif
            enddo

            pointRight(:,i,(j-1)*(interpNum+1)+1) = pointLeft(:,i,(number-1)*(interpNum+1)+1)
        
        enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  Leftの補間
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,nodeNum

        !!中心点を計算
        !
        ! 本当は全nodeからの距離と距離の標準偏差が最小になる位置を探す方が良い
        !
        center(:) = 0.0d0
  
        do j = 1, inputWNum
            center(:) = center(:) + pointLeft(:, i, (j - 1)*(interpNum + 1) + 1)
        enddo
        center(:) = center(:) / dble(inputWNum)

        !!!!!  wireを中心に向かって移動（New）  !!!!!
        do j = 1, inputWNum
            adjst(:) = center(:) - pointLeft(:, i, (j - 1)*(interpNum + 1) + 1)
            temp = sqrt(adjst(1)**2 + adjst(2)**2 + adjst(3)**2)
            adjst(:) = adjst(:) / temp
            pointLeft(:,i,(j-1)*(interpNum + 1)+1) = pointLeft(:,i,(j-1)*(interpNum + 1)+1) + ds * adjst(:)
        enddo

        !!!!!  wireを補間  !!!!!
        do j = 1,inputWNum
            p(:) = pointLeft(:,i,(j-1)*(interpNum+1)+1) - center(:)
            if(j .ne. wireNum/(interpNum+1))then
                q(:) = pointLeft(:,i,(j)*(interpNum+1)+1) - center(:)
            else
                q(:) = pointLeft(:,i,1) - center(:)
            endif

            !!円の半径を計算
            radii1 = sqrt(p(1)**2 + p(2)**2 + p(3)**2 )
            radii2 = sqrt(q(1)**2 + q(2)**2 + q(3)**2 )

            !!外積から中心点を通る法線を計算
            !  ＜将来的な候補＞楕円を仮定して、最小二乗法より式を得る方法
  
            !!外積計算  point((j-1)*4+1) × point((j)*4+1)

            Normal(1) = p(2) * q(3) - p(3) * q(2)
            Normal(2) = p(3) * q(1) - p(1) * q(3)
            Normal(3) = p(1) * q(2) - p(2) * q(1)
            temp = sqrt(Normal(1)**2 + Normal(2)**2 + Normal(3)**2)
            Normal(1) = Normal(1)/temp
            Normal(2) = Normal(2)/temp
            Normal(3) = Normal(3)/temp

            theta = acos( dot_product(p,q) / (radii1 * radii2) )

            !!点を回転
            l(:) = Normal(:)

            do k = 1,interpNum !!補間
                pointLeft(1,i,(j-1)*(interpNum+1)+1+k) = (l(1)**2 + (1.0d0 - l(1)**2)*cos(theta/4.0d0*dble(k))          ) * p(1) + &
                           (l(1)*l(2) * (1.0d0 - cos(theta/4.0d0*dble(k)))-l(3)*sin(theta/4.0d0*dble(k))) * p(2) + &
                           (l(1)*l(3) * (1.0d0 - cos(theta/4.0d0*dble(k)))+l(2)*sin(theta/4.0d0*dble(k))) * p(3)
                           
                pointLeft(2,i,(j-1)*(interpNum+1)+1+k) = (l(1)*l(2) * (1.0d0 - cos(theta/4.0d0*dble(k))) + l(3) * &
                                                                          sin(theta/4.0d0*dble(k)))   * p(1) + &
                           (l(2)**2 + (1.0d0 - l(2)**2)*cos(theta/4.0d0*dble(k)))                 * p(2) + &
                           (l(2)*l(3) * (1.0d0 - cos(theta/4.0d0*dble(k)))-l(1)*sin(theta/4.0d0*dble(k))) * p(3)
                           
                pointLeft(3,i,(j-1)*(interpNum+1)+1+k) = (l(1)*l(3) * (1.0d0 - cos(theta/4.0d0*dble(k))) - l(2) * &
                                                                          sin(theta/4.0d0*dble(k)))   * p(1) + &
                           (l(2)*l(3) * (1.0d0 - cos(theta/4.0d0*dble(k)))+l(1)*sin(theta/4.0d0*dble(k))) * p(2) + &
                           (l(3)**2 + (1.0d0 - l(3)**2)*cos(theta/4.0d0*dble(k)))                 * p(3)
                           
                !!半径を調節
                rnew = dble(k) / dble(interpNum+1) * (radii2 - radii1) + radii1
                pointLeft(:,i,(j-1)*(interpNum+1)+1+k) = pointLeft(:,i,(j-1)*(interpNum+1)+1+k) / radii1 * rnew + center(:)

            enddo
        enddo
    enddo










end program
