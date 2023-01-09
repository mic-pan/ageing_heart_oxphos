using VegaLite

function scale_svg(svg,scale=1.5)
    wstr = "width"
    hstr = "height"
    
    (iw,width) = find_length(svg,wstr)
    (ih,height) = find_length(svg,hstr)
    
    new_width = width*scale
    new_height = height*scale
    
    wrule = ("$wstr=\"$(svg[iw])\"" => "$wstr=\"$new_width\"")
    hrule = ("$hstr=\"$(svg[ih])\"" => "$hstr=\"$new_height\"")
    
    svg1 = replace(svg,wrule;count=1)
    svg2 = replace(svg1,hrule;count=1)
end

function enlarge(v::VegaLite.VLSpec, scale=1.5)
    svg = scale_svg(VegaLite.convert_vl_to_svg(v),scale)
    display(MIME("image/svg+xml"), svg)
end

function find_length(svg,str)
    i1 = findfirst(str*"=\"",svg)[end]+1
    i2 = findnext("\"",svg,i1+1)[end]-1
    ilen = i1:i2
    len = parse(Float64,svg[ilen])
    return ilen,len
end

function savevl(fig,name,formats=["svg","pdf"])
    for f in formats
        fig |> save(name*"."*f)
    end
end