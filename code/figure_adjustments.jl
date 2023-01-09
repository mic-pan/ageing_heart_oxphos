function set_axis_title!(fig,axis,title)
    fig.layer[1]["encoding"][axis]["axis"]["title"] = title
end

function set_domain!(fig,axis,lims)
    fig.layer[1]["encoding"][axis]["scale"]["domain"] = lims
    mark = fig.layer[1]["mark"]
    if mark isa AbstractDict
        mark["clip"] = true
    else
        fig.layer[1]["mark"] = OrderedDict(
            "type" => mark,
            "clip" => true
        )
    end
end

function zero_axis!(fig,axis,zero)
    fig.layer[1]["encoding"][axis]["scale"]["zero"] = zero
end
