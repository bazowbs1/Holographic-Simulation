function OutPicture = CustomizedZoom( InPicture,ZoomFactor )

    ySize = size( InPicture,1 ); yCrop = floor( ySize*abs(ZoomFactor-1 )/2 );
    xSize = size( InPicture,2 ); xCrop = floor( xSize*abs(ZoomFactor-1 )/2 );
    zSize = size( InPicture,3 );
    
    if ZoomFactor == 0
        OutPicture = double( zeros( ySize,xSize,zSize ) );
        return;
    end
    
    if ZoomFactor == 1
        OutPicture = InPicture;
        return;
    end
        
    ZoomPicture = imresize( InPicture,ZoomFactor,'cubic' );
    yySize = size( ZoomPicture,1 );
    xxSize = size( ZoomPicture,2 );
    
    if ZoomFactor>1 %shrink
        OutPicture = ZoomPicture( yCrop+1:yCrop+ySize,xCrop+1:xCrop+xSize,: );
    else %expand
        OutPicture = double(zeros(ySize,xSize,zSize));
        OutPicture( yCrop+1:yCrop+yySize,xCrop+1:xCrop+xxSize,: ) = ZoomPicture;
    end
    
end
