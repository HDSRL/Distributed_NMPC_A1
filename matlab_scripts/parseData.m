function [varargout] = parseData(filenames)
    tailIncluded = true;
    if ~iscell(filenames)
        filenames = {filenames};
        tailIncluded = false;
    elseif length(filenames)==1
        tailIncluded = false;
    end
    
    gearReduction = repmat([1,1,1],1,4);

    V60.raw = csvread(filenames{1});
    V60.time = V60.raw(:,1);
    V60.outputs = V60.raw(:,2:13);
    V60.torques = V60.raw(:,14:25)./gearReduction;
    V60.pos = V60.raw(:,26:43);
    V60.vel = V60.raw(:,44:61);
    V60.hd = V60.raw(:,62:73);
    V60.dhd = V60.raw(:,74:85);
    V60.ddhd = V60.raw(:,86:97);
    V60.trajMPC = V60.raw(:,98:109);
    V60.forceMPC = V60.raw(:,110:121);
    V60.doutputs = V60.raw(:,122:133);
    V60.V = V60.raw(:,134);
    V60.dV = V60.raw(:,135);
    V60.force = V60.raw(:,136:139);
    V60.phase = V60.raw(:,140);
%     V60.QPforce = V60.raw(:,140:151);
    
    varargout{1} = V60;

    if tailIncluded
        tail.raw = csvread(filenames{2});
        tail.time = tail.raw(:,1);
        tail.domain = tail.raw(:,2);
        tail.outputs = tail.raw(:,3:4);
        tail.torques = tail.raw(:,5:6);
        tail.indivTorque = tail.raw(:,7:15);
        tail.pos = tail.raw(:,16:17);
        tail.vel = tail.raw(:,18:19);
        varargout{2} = tail;
    end
end

