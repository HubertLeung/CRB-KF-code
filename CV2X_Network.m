function  [posAP,sign,numAP] = CV2X_Network(nodeLocUpdate,dim,squSize,delPosAP)

networkLim = zeros(dim,2);
for idx = 1:dim
    networkLim(idx,1) = min(nodeLocUpdate(dim+idx,:));
    networkLim(idx,2) = max(nodeLocUpdate(dim+idx,:));
end


for idx = 1:dim
    networkLim(idx,1) = floor(networkLim(idx,1))-delPosAP;
end

networkLim(1,2) = networkLim(1,1) + squSize;
networkLim(2,2) = networkLim(2,1) + squSize;
networkLim(3,2) = networkLim(3,1) + 2*delPosAP;

dim1 = (networkLim(1,1):delPosAP:networkLim(1,2));
dim2 = (networkLim(2,1):delPosAP:networkLim(2,2));
dim3 = (networkLim(3,1):delPosAP:networkLim(3,2));

L1 = length(dim1);
L2 = length(dim2);
L3 = length(dim3);

numAP = L1*L2*L3;
posAP = zeros(dim,numAP);

sign = 1;
for idx1 = 1:L1
    for idx2 = 1:L2
        for idx3 = 1:L3
            posAP(:,sign) = [dim1(idx1) dim2(idx2) dim3(idx3)]';
            sign = sign + 1;
        end
    end
end