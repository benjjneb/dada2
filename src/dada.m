% DADA, the Divisive Amplicon Denoising Algorithm
%     Copyright (C) 2012 Michael Jeremy Rosen
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function dada(inFolder,varargin)
p = inputParser;
p.addRequired('inFolder',@ischar);
p.addOptional('omegaA',.01,@(x)isnumeric(x) && isscalar(x) && x>=0);
p.addOptional('omegaR',.01,@(x)isnumeric(x) && isscalar(x) && x>=0);
%G: standard gap-open and gap-extension penalty
p.addParamValue('G',-4,@(x)isnumeric(x) && isscalar(x));
%GH: homopolymer gap penalty
p.addParamValue('GH',-1,@(x)isnumeric(x) && isscalar(x));
%low: the threshold for the largest sequence in a cluster above which it be
%taken to be the consensus sequence of the cluster. if no sequence in a
%cluster exceeds this value, a multiple alignment will be used to generate
%a consensus sequence
p.addParamValue('low',3,@(x)isnumeric(x) && isscalar(x));
%context: T/F value determining whether to use nearest-neighbor
%context-dependent error probabilities when computing lambdas
p.addParamValue('context',false,@(x)islogical(x) && isscalar(x));
%err: the filename of an err matrix file
p.addParamValue('err','',@ischar);
%noUpdate: run DADA with no error rate updating. useful in combination with
%an imported error rate matrix already believed to be correct
p.addParamValue('noUpdate',false,@(x)islogical(x) && isscalar(x));
%verbose: T/F value determining whether to run in verbose mode
p.addParamValue('verbose',false,@(x)islogical(x) && isscalar(x));
%tol: tolerance for norms in T convergence
p.addParamValue('tol',1e-9,@(x)isnumeric(x) && isscale(x));
p.parse(inFolder,varargin{:});
inFolder = p.Results.inFolder;
omegaA = p.Results.omegaA;
omegaR = p.Results.omegaR;
G = p.Results.G;
GH = p.Results.GH;
LOW_ABUNDANCE_THRESHOLD = p.Results.low;
context = p.Results.context;
err = p.Results.err;
noUpdate = p.Results.noUpdate;
VERBOSE = p.Results.verbose;
tol = p.Results.tol;
if G ~= GH
    HOMO = true;
else
    HOMO = false;
end
prevERR = zeros(100,16); %matrix of average error rates on previous rounds
%create a name for the output folder
if ~HOMO
    outFolder = [datestr(now,30) '_in=' inFolder '_omegaA=' ...
        num2str(omegaA) '_omegaR=' num2str(omegaR) '_context = ' ...
        num2str(context) '_G=' num2str(G)];
else
    outFolder = [datestr(now,30) '_in=' inFolder '_omegaA=' ...
        num2str(omegaA) '_omegaR=' num2str(omegaR) '_context = ' ...
        num2str(context) '_G=' num2str(G) '_GH=' num2str(GH)];
end
mkdir(outFolder);
amps = dir(inFolder); %the amplicons
amps = amps([amps.bytes]~=0); %remove empty files

if ~isempty(err) %user passed in error probabilities
    iteration = 1;
    load(err,'ERRa','ERR');
else
    iteration = 0;
    %the 0th iteration does not include any clustering: each input file
    %will result in a single output bin. the only purpose is to estimate
    %initial error probabilities.
    ERRa = ones(4,4);
    ERR = ones(4,4,4,4);
    scoreMat = nuc44;
    scoreMat = scoreMat(1:4,1:4);
end

%termination conditions:
%1. we have already clustered with the current value of T
%2. we were told never to update T after the first iteration
%3. the norm of the difference of T on consecutive rounds falls below tol

while ~ismember(ERRa(:)',prevERR,'rows') && (~noUpdate || iteration==1) ...
        && ~( iteration > 1 && ...
        norm(ERRa(:)'-prevERR(iteration-1,:)) < tol)
    
    if iteration > 10
        disp 'Run for 10 iterations. Quitting.';
        break;
    end
    
    mkdir(outFolder,num2str(iteration)); %new output folder
    if iteration > 0
        prevERR(iteration,:) = ERRa(:);
        scoreMat = 5 + log(ERRa);
    end
    disp 'current score matrix:';
    disp(scoreMat);
    if min(scoreMat(:)) < 2 * G
        disp(['BEWARE: some elements of the score matrix are less than'...
            ' twice the gap penalty, ' num2str(G) '. All substitutions'...
            ' of this type will therefore be completely gapped out and'...
            ' the error rate will collapse to 0. Consider exiting and'...
            ' boosting the size of the gap penalty.']);
    end
    %loop clustering the sequences in each input file starts here
    for amp = 1:length(amps)
        disp(amps(amp).name);
        fid = fopen([inFolder '/' amps(amp).name],'r');
        tmp = textscan(fid,'%f\t%s');
        fclose(fid);
        reads = tmp{1};
        seqs = tmp{2};
        %size order input seqs        
        [reads,I] = sort(reads,'descend');
        seqs = seqs(I);

        %remove seqs with 'N's (ambig calls)
        for seq = length(seqs):-1:1
            if ~isempty(strfind(seqs{seq},'N')) || ...
                    ~isempty(strfind(seqs{seq},'n'))
                seqs(seq) = [];
                reads(seq) = [];
                continue;
            end
        end
        %convert sequences to uint8 format, used throughout the algorithm
        for seq = 1:length(seqs)
            seqs{seq} = uint8(nt2int(seqs{seq}));
        end        
        if amp == 1
            %initialize struct for singleton pval computations. the
            %60/50/40 suffixes on the degeneracy, dot product, and cumsum
            %fields represent ideal templates with different GC contexts.
            %currently, no preparations are made to handle templates with
            %different lengths. with no computational constrains, psingle
            %could include degeneracy vectors for each cluster.
            E = struct('dmax',[],'L',[],'e1',[],'e2',[],'multi',[],...
                'gc',struct('gcfrac',[],'bc',[],'bin',[],'self',[],'m',[],...
                'dot',[],'p',[]),'Lmaxdmax',[]);
            E.e1 = diag(ERRa); %4x1 of non-error probs
            %4x3 of error probs: what must multiply self by to get rate
            E.e2 = [[ERRa(1,2) ERRa(1,3) ERRa(1,4)] / ERRa(1,1); ...
                [ERRa(2,1) ERRa(2,3) ERRa(2,4)] / ERRa(2,2); ...
                [ERRa(3,1) ERRa(3,2) ERRa(3,4)] / ERRa(3,3); ...
                [ERRa(4,1) ERRa(4,2) ERRa(4,3)] / ERRa(4,4)];
            %initialize E for dmax = 0
            E.dmax = 0;
            %E.L is the collection of ordered rates out to distance dmax
            %including the no-error sequence. These are relative error
            %rates, meaning that they give the true error rate of the
            %sequence when multiplied by the selfing rate. e.g. the E.L
            %entry for being exactly right has value 1. Relative rates are
            %used so that different GC-contents, which have different
            %selfing rates, can share a single set of rates against which
            %query lambda are compared
            E.L = 1;
            %E.Lmaxdmax is the largest error rate seen at distance dmax. if
            %a query lambda arrives that is smaller than this, than there
            %is a possibility that there are larger lambda a larger
            %distance, so this is the condition for triggering expanding
            %dmax
            E.Lmaxdmax = 1;
            %the length of the most abundant sequence in first cluster is
            %assumed for p-value computations for now
            ampL = length(seqs{1}); %assumed length during p-value comp
            gcfrac = [.6 .5 .4]; %GC-fraction to interpolate between
            for ii = 1:length(gcfrac)
                %compute a basecount consistent with each GC-frac and ampL
                gc = ceil(ampL*gcfrac(ii)); g = ceil(gc/2); c = gc - g;
                at = length(seqs{1}) - gc; a = ceil(at/2); t = at - a;
                E.gc(ii).gcfrac = gcfrac(ii);
                E.gc(ii).bc = [a;c;g;t];
                self = prod(E.e1.^[a;c;g;t]); %selfing rate
                E.gc(ii).self = self;
                E.gc(ii).m = 1;
                E.gc(ii).dot = self;
                E.gc(ii).p = self;
            end
        end
        
        %bin: the struct containing all information about each cluster
        bin = struct('seq',[],'R',[],'baseCount',[],'self',[],...
            'famHash',java.util.Hashtable,'fam',struct('seq',[],'r',[],...
            'lambda',[],'p',[],'key','','raw',struct('seq',{},...
            'reads',[],'lambda',[],'e',[],'subPos',{},'subNT',{},...
            'subkey',{})),'pS',[],'pSf',[],'pSval',[]);
        elsUpdated = []; %flag to recompute R & template
        RUpdated = []; %flag to recompute pvals
        templateUpdated = []; %flag to recompute lambdas & pvals
        budB = []; %global variable set by sig and used by bud
        budF = []; %global variable set by sig and used by bud
        
        %innermost loop, adding clusters until reaching statistical
        %consistency, beings here:
        init(); %place all sequences in bin 1
        templateUpdate();
        RUpdate();
        lambdaUpdate();
        famUpdate();
        if iteration > 0
            eUpdate();
            pUpdate();
            %iterate until no sequence changes bins due to abundance
            while sig();
                bud();
                while ~isempty(elsUpdated) %some els moved
                    templateUpdate();
                    RUpdate();
                    lambdaUpdate();
                    famUpdate();
                    eUpdate();
                    pUpdate();
                    shuffle();
                end
            end
        end
        %dump output to file
        if regexp(amps(amp).name,'[.]') %filenames has extension
            fname = [outFolder '/' num2str(iteration) '/' ...
                regexprep(amps(amp).name,'[^.]+$','mat')]; %replace w .mat
        else %no extension
            fname = [outFolder '/' num2str(iteration) '/' ...
                amps(amp).name '.mat']; %just append .mat
        end
        out(bin,fname,scoreMat,LOW_ABUNDANCE_THRESHOLD);
    end
    [ERRa,ERR] = MLE(outFolder,iteration); %reestimate error rates
    iteration = iteration + 1;
end
    function init()
        %seeds clustering by placing all sequences in one cluster, bin(1)
        bin(1).R = 0; %this is set to 0 so that N_updated is triggered
        [R,I] = max(reads);
        if R <= LOW_ABUNDANCE_THRESHOLD
            bin(1).seq = consensus(scoreMat,seqs,reads);
        else
            bin(1).seq = seqs{I};
        end
        %place all seqs in family 1
        for i = 1:length(seqs)
            bin(1).fam(1).raw(i).seq = seqs{i};
            bin(1).fam(1).raw(i).reads = reads(i);
        end
        elsUpdated = 1;
    end

    function tf = sig()
        %determines whether statistical consistency has been reached
        tf = false;
        %find bin and fam of minimum abundance p-value
        minp = 1;
        for b = 1:length(bin)
            [p,f] = min([bin(b).fam.p]);
            if p < minp
                B = b;
                F = f;
                minp = p;
            end
        end
        if minp * length([bin.fam]) < omegaA
            budB = B;
            budF = F;
            tf = true;
        else %omegaR values are stored via Boolean values.
            lambdaMin = inf;
            for b = find([bin.pS]);
                if bin(b).fam(bin(b).pSf).lambda < lambdaMin
                    lambdaMin = bin(b).fam(bin(b).pSf).lambda;
                    budB = b;
                    budF = bin(budB).pSf;
                    tf = true;
                end
            end
        end
    end

    function bud()
        %starts a new cluster
        B = budB;
        F = budF;
        if VERBOSE
            disp(['starting new bin with family ' num2str(F) ...
                ' from bin ' num2str(B) ...
                ',lambda:' num2str(bin(B).fam(F).lambda) ...
                ',R:' num2str(bin(B).R) ...
                ',p:' num2str(bin(B).fam(F).p) ...
                ',self:' num2str(bin(B).self)]);
        end
        %remove the fam from its old bin and seed new one
        bin(end+1).fam(1) = bin(B).fam(F);
        bin(B).fam(F) = [];
        
        %remove from old hash table and fix remaining values
        bin(B).famHash.remove( bin(end).fam(1).key );
        for f = F:length(bin(B).fam)
            bin(B).famHash.put( bin(B).fam(f).key , f );
        end
        
        [R,I] = max([bin(end).fam(1).raw.reads]);
        if R <= LOW_ABUNDANCE_THRESHOLD
            bin(end).seq = consensus(scoreMat,{bin(end).fam(1).raw.seq},...
                [bin(end).fam(1).raw.reads]);
        else
            bin(end).seq = bin(end).fam(1).raw(I).seq;
        end
        bin(end).R = 0; %this is set to 0 so that updates are triggered
        elsUpdated = [B length(bin)];
    end

    function templateUpdate()
        %updates the consensus sequence
        templateUpdated = zeros(1,length(bin));
        for b = elsUpdated
            if bin(b).R > 0 %not a new bin
                maxr = 0;
                for f = 1:length(bin(b).fam)
                    [r,ind] = max([bin(b).fam(f).raw.reads]);
                    if r > maxr
                        maxr = r;
                        F = f;
                        R = ind;
                    end
                end
                if maxr <= LOW_ABUNDANCE_THRESHOLD
                    raws = [bin(b).fam.raw]; %concatenate raws together
                    t = consensus(scoreMat,{raws.seq},[raws.reads]);
                else
                    t = bin(b).fam(F).raw(R).seq;
                end
            else %new bin. don't reupdate the template
                t = bin(b).seq;
            end
            if ~isequal(t,bin(b).seq) || bin(b).R == 0
                %update seq and baseCount fields
                bin(b).seq = t;
                for nt = 1:4
                    bin(b).baseCount(nt) = sum(t==nt);
                end
                %update the self-production lambda.
                t = t ( t ~= 5 );
                lambda = 1;
                for i = 2:length(t)-1
                    if context
                        lambda = lambda * ERR(t(i-1),t(i),t(i),t(i+1));
                    else
                        lambda = lambda * ERRa(t(i),t(i));
                    end
                end
                bin(b).self = lambda;
                if VERBOSE
                    disp(['replacing bin ' num2str(b) ' template with']);
                    disp(int2nt(bin(b).seq));
                    disp(['and self-lambda ' num2str(lambda)]);
                end
                templateUpdated(b) = 1;
            end
        end
        templateUpdated = find(templateUpdated);
    end

    function RUpdate()
        %updates R, the total number of reads in the bin
        RUpdated = zeros(1,length(bin));
        for b = elsUpdated
            R = 0;
            for f = 1:length(bin(b).fam)
                for r = 1:length(bin(b).fam(f).raw)
                    R = R + bin(b).fam(f).raw(r).reads;
                end
            end
            if  R ~= bin(b).R
                bin(b).R = R;
                RUpdated(b) = 1;
            end
        end
        RUpdated = find(RUpdated);
    end

    function lambdaUpdate()
        %update the alignments and lambda of all raws to new consensus
        %sequences
        for i = templateUpdated
            for j = 1:length(bin)
                for f = 1:length(bin(j).fam)
                    for r = 1:length(bin(j).fam(f).raw)
                        if ~HOMO
                            [subPos, subNT] = align2dom( bin(i).seq, ...
                                bin(j).fam(f).raw(r).seq, scoreMat, G);
                            bin(j).fam(f).raw(r).subPos{i} = subPos;
                            bin(j).fam(f).raw(r).subNT{i} = subNT;
                            bin(j).fam(f).raw(r).subkey{i} = ...
                                [sprintf('%d,',subPos) sprintf('%d,',subNT)];
                        else
                            [subPos, subNT] = align2dom_homo( bin(i).seq, ...
                                bin(j).fam(f).raw(r).seq, scoreMat, G, GH);
                            bin(j).fam(f).raw(r).subPos{i} = subPos;
                            bin(j).fam(f).raw(r).subNT{i} = subNT;
                            bin(j).fam(f).raw(r).subkey{i} = ...
                                [sprintf('%d,',subPos) sprintf('%d,',subNT)];
                        end
                        bin(j).fam(f).raw(r).lambda(i) = lambdaShift( ...
                            bin(i).seq, bin(j).fam(f).raw(r).subPos{i}, ...
                            bin(j).fam(f).raw(r).subNT{i}, ...
                            bin(i).self, ERR, ERRa, context );
                        if VERBOSE
                            disp(['b:' num2str(j) ' f:' num2str(f) ' r:'...
                                num2str(r) ' lambda(' num2str(i) '):'...
                                num2str(bin(j).fam(f).raw(r).lambda(i)) ...
                                ' subkey{' num2str(i) '}:' ...
                                bin(j).fam(f).raw(r).subkey{i} ' reads:'...
                                num2str(bin(j).fam(f).raw(r).reads)]);
                            
                        end
                    end
                end
            end
        end
    end

    function famUpdate()
        %reform families for bins with updated templates
        for i = templateUpdated
            raws = [bin(i).fam.raw]; %concatenate all the raws together
            bin(i).famHash = java.util.Hashtable;
            bin(i).fam = [];
            for j = 1:length(raws)
                f = bin(i).famHash.get( raws(j).subkey{i} );
                if ~isempty(f) %join existing family, f
                    bin(i).fam(f).raw(end+1) = raws(j);
                    bin(i).fam(f).r = bin(i).fam(f).r + raws(j).reads;
                else
                    f = length(bin(i).fam) + 1; %start new family, f
                    bin(i).fam(f).raw(1) = raws(j);
                    bin(i).famHash.put( raws(j).subkey{i}, f);
                    bin(i).fam(f).key = raws(j).subkey{i};
                    bin(i).fam(f).seq = bin(i).seq;
                    bin(i).fam(f).seq(raws(j).subPos{i}) = ...
                        raws(j).subNT{i};
                    bin(i).fam(f).r = raws(j).reads;
                    bin(i).fam(f).lambda = raws(j).lambda(i);
                end
            end
            if VERBOSE
                disp(['reformed ' num2str(length(bin(i).fam)) ...
                    ' families for bin ' num2str(i) ...
                    ' with read numbers:']);
                disp([bin(i).fam.r]);
            end
        end
    end

    function eUpdate()
        %update expected number of reads of each sequence
        for i = union(templateUpdated,RUpdated)
            for j = 1:length(bin)
                for f = 1:length(bin(j).fam)
                    for r = 1:length(bin(j).fam(f).raw)
                        bin(j).fam(f).raw(r).e(i) = ...
                            bin(j).fam(f).raw(r).lambda(i) * bin(i).R;
                        if VERBOSE
                            disp(['b:' num2str(j) ' f:' num2str(f) ' r:'...
                                num2str(r) ' E_reads(' num2str(i) '):'...
                                num2str(bin(j).fam(f).raw(r).e(i))]);
                        end
                    end
                end
            end
        end
    end

    function pUpdate()
        for b = union(templateUpdated,union(RUpdated,elsUpdated))
            %update the abundance p-values.
            R = bin(b).R;
            for f = 1:length(bin(b).fam)
                r = bin(b).fam(f).r;
                lambda = bin(b).fam(f).lambda;
                if isequal(bin(b).seq,bin(b).fam(f).seq)
                    bin(b).fam(f).p = 1; %bin center
                elseif r == 1
                    bin(b).fam(f).p = 1; %singleton
                else
                    %compute normalization
                    norm = 1 - poisspdf(0,R*lambda);
                    if norm == 0
                        norm = lambda*R; %better approximation
                    end
                    tail = 1 - poisscdf(r-1,R*lambda);
                    if tail == 0
                        tail = poisspdf(r,R*lambda); %use first term
                    end
                    bin(b).fam(f).p = tail / norm;
                end
            end
            Nfam = 0;
            for i = 1:length(bin)
                Nfam = Nfam + length(bin(i).fam);
            end
            if min([bin(b).fam.p]) * Nfam < omegaA
                continue; %there is a significant abundance pval
            end
            %update the singleton p-value.
            [lambda,bin(b).pSf] = min([bin(b).fam.lambda]/bin(b).self);
            if lambda == 1 %there is only one family.
                bin(b).pS = false;
                continue;
            end
                        
            while(1)
                %compute a pval
                x = find(E.L > lambda,1,'last');
                fitX = zeros(1,length(E.gc)); %the GC content
                fitY = zeros(1,length(E.gc)); %the p-value
                for i = 1:length(E.gc)
                    fitX(i) = E.gc(i).gcfrac;
                    fitY(i) = log10(E.gc(i).p(x));
                end
                %perform least-squares linear fit
                fit = polyfit(fitX,fitY,1);
                %interpolate the log10 of the p-value for the observed 
                %lambda and transform back
                p =  bin(b).R * 10^(fit(2) + fit(1) * ...
                    (sum(bin(b).baseCount([2 3])) / length(bin(b).seq)));
%                 p = 1 - (1-10^(fit(2)+fit(1)*...
%                     (sum(bin(b).baseCount([2 3]))/...
%                     length(bin(b).seq))))^bin(b).R;
                bin(b).pSval = p;
                if p * length(bin) < omegaR
                    %                 if p * sum(arrayfun(@(s)length(s.fam),bin) > 1) < omegaR
                    %lambda is significant even without expanding dmax
                    bin(b).pS = true;                    
                    break;
                elseif lambda < E.Lmaxdmax && E.dmax < 10  % expanding pval vectors
                    %expand dmax and try computing p again
                    d = E.dmax + 1;
                    disp(['expanding E to distance ' num2str(d)]);
                    E.dmax = d;
                    %expand multi cell array. this is used to lookup
                    %multinomial coefficients when computing error degeneracies
                    E.multi{d} = zeros(d+1,d+1,d+1);
                    for i = ceil(d/3):d
                        for j = ceil((d-i)/2):min(i,d-i)
                            k = d - i - j;
                            n = factorial(d) / prod(factorial([i j k]));
                            E.multi{d}(i+1,j+1,k+1) = n;
                            E.multi{d}(i+1,k+1,j+1) = n;
                            E.multi{d}(j+1,i+1,k+1) = n;
                            E.multi{d}(j+1,k+1,i+1) = n;
                            E.multi{d}(k+1,i+1,j+1) = n;
                            E.multi{d}(k+1,j+1,i+1) = n;
                        end
                    end
                    %expand bin matrices. these are used to lookup binomial
                    %coefficients when computing error degeneracies.
                    %E.gc(i).bin(d,j) is number of ways of putting d errors
                    %onto the bases of type j of GC-content i.
                    for i = 1:length(E.gc)
                        for j = 1:4
                            if d <= E.gc(i).bc(j)
                                E.gc(i).bin(d,j) = nchoosek(E.gc(i).bc(j),d);
                            else
                                E.gc(i).bin(d,j) = 0;
                            end
                        end
                    end
                    %compute all error probabilities and degeneracies for each
                    %gc-fraction at new distance d and store in L vector and m
                    %cell array of vectors of same length as L.
                    %nchoosek(d+11,11) is the number of distinct error types at
                    %distance d.
                    L = zeros(nchoosek(d+11,11),1);
                    m = cell(1,length(E.gc));
                    for i = 1:length(m)
                        m{i} = ones(nchoosek(d+11,11),1);
                    end
                    %Let A be a matrix that keeps track of each error type:
                    %A(i,j) contains the number of i->j errors. we are
                    %considering a distance d so the elements of A always
                    %add to d. A simple procedure is used within the while loop
                    %that makes sure that we iterate through all matrices whose
                    %elements add to d. we keep an index value i throughout.
                    A = zeros(4,3);
                    %note that linear indexing of A is used throughout
                    A(1) = d; %start with d A->C errors
                    i = 1;
                    while(1)
                        bases = sum(A,2); %bases involved in errors
                        L(i) = prod(prod(E.e2.^A)); %relative error rate of A
                        for j = 1:length(m) %compute the degeneracy of A
                            for k = 1:4
                                if bases(k) > 0
                                    %multiply number of sets of sites that
                                    %could be occupied by errors away from k
                                    m{j}(i) = m{j}(i) * E.gc(j).bin(bases(k),k);
                                    %multiply number of distinct ways to put
                                    %errs away from k onto these sites -- i.e.
                                    %ways to rearrange the three different
                                    %error types away from this base
                                    m{j}(i) = m{j}(i) * E.multi{bases(k)}(...
                                        A(k,1)+1,A(k,2)+1,A(k,3)+1);
                                end
                            end
                        end
                        if A(end) == d %d T->G errors is termination condition
                            break;
                        elseif A(end) == 0 %no T->G errs: advance final error
                            x = find(A,1,'last');
                            A(x) = A(x) - 1;
                            A(x+1) = A(x+1) + 1;
                        else
                            %some T->G errs. move them back to closest error,
                            %pick up on error from that error, and advance
                            %these all forward 1
                            n = A(end);
                            A(end) = 0;
                            x = find(A,1,'last');
                            A(x) = A(x) - 1;
                            A(x+1) = n + 1;
                        end
                        i = i + 1; %advance index
                    end 
                    E.Lmaxdmax = max(L);
                    %interpolate L into existing error rates
                    [E.L,I]=sort([E.L;L],'descend');
                    for i = 1:length(E.gc)
                        E.gc(i).m = [E.gc(i).m;m{i}];
                        E.gc(i).m = E.gc(i).m(I);
                        E.gc(i).dot = [E.gc(i).dot;E.gc(i).self*L.*m{i}];
                        E.gc(i).dot = E.gc(i).dot(I);
                        E.gc(i).p = 1 - cumsum(E.gc(i).dot);
                    end
                else
                    %p won't be changed by expanding dmax and lambda is not
                    %statistically significant
                    bin(b).pS = false;
                    break;
                end
            end
        end
    end

    function shuffle()
        %move each sequence to the bin that produces the highest expected
        %number of that sequence
        elsUpdated = [];
        for b = 1:length(bin)
            fam_remove = false(1,length(bin(b).fam));
            for f = 1:length(bin(b).fam)
                if isequal(bin(b).seq,bin(b).fam(f).seq)
                    continue;
                end
                raw_remove = false(1,length(bin(b).fam(f).raw));
                for r = 1:length(bin(b).fam(f).raw)
                    [~,c] = max(bin(b).fam(f).raw(r).e);
                    if c ~= b %raw should move from bin 'b' to bin 'c'
                        if VERBOSE
                            disp(['moving raw ' num2str(r) ...
                                ' from family ' num2str(f) ...
                                ' with ' ...
                                num2str(bin(b).fam(f).raw(r).reads) ...
                                ' reads of bin ' num2str(b) ...
                                '(where it had lambda: ' ...
                                num2str(bin(b).fam(f).raw(r).lambda(b)) ...
                                ' e: ' num2str(bin(b).fam(f).raw(r).e(b))...
                                ') to bin ' num2str(c) ...
                                '(where it will have lambda: ' ...
                                num2str(bin(b).fam(f).raw(r).lambda(c)) ...
                                ' e: ' num2str(bin(b).fam(f).raw(r).e(c))...
                                ')']);
                        end
                        bin(b).fam(f).r = bin(b).fam(f).r - ...
                            bin(b).fam(f).raw(r).reads;
                        %check whether to join family or start new one
                        F = bin(c).famHash.get( ...
                            bin(b).fam(f).raw(r).subkey{c} );
                        if ~isempty(F)
                            if VERBOSE
                                disp(['joining family ' num2str(F)]);
                            end
                            bin(c).fam(F).raw(end+1) = ...
                                bin(b).fam(f).raw(r);
                            bin(c).fam(F).r = bin(c).fam(F).r + ...
                                bin(b).fam(f).raw(r).reads;
                        else
                            F = length(bin(c).fam) + 1;
                            if VERBOSE
                                disp(['starting new family ' num2str(F)]);
                            end
                            bin(c).fam(F).raw(1) = bin(b).fam(f).raw(r);
                            bin(c).famHash.put(...
                                bin(b).fam(f).raw(r).subkey{c}, F);
                            bin(c).fam(F).key = ...
                                bin(b).fam(f).raw(r).subkey{c};
                            bin(c).fam(F).seq = bin(c).seq;
                            bin(c).fam(F).seq( ...
                                bin(b).fam(f).raw(r).subPos{c} ) = ...
                                bin(b).fam(f).raw(r).subNT{c};
                            bin(c).fam(F).r = bin(b).fam(f).raw(r).reads;
                            bin(c).fam(F).lambda = ...
                                bin(b).fam(f).raw(r).lambda(c);
                        end
                        raw_remove(r) = true;
                        if isempty(elsUpdated)
                            %this extra condition is needed because as of
                            %R2013a, the union of a row vector is an empty
                            %vector is a column vector.
                            elsUpdated = [b c];
                        else
                            elsUpdated = union(elsUpdated,[b c]);
                        end
                    end
                end
                bin(b).fam(f).raw(raw_remove) = [];
                if isempty(bin(b).fam(f).raw)
                    fam_remove(f) = true;
                end
            end
            %remove keys of empties from hash table
            for f = find(fam_remove)
                bin(b).famHash.remove( bin(b).fam(f).key );
            end
            %remove empty families
            bin(b).fam(fam_remove) = [];
            %update hashtable values of non-empty families
            bin(b).famHash = java.util.Hashtable;
            for f = 1:length(bin(b).fam)
                bin(b).famHash.put( bin(b).fam(f).key, f);
            end
        end
    end
end

function out(bin,fname,scoreMat,LOW_ABUNDANCE_THRESHOLD)
%dump the bin struct array and some other information to a .mat file
reals = cell(1,length(bin));
for i = 1:length(bin)
    maxr = 0;
    for j = 1:length(bin(i).fam)
        [r,ind] = max([bin(i).fam(j).raw.reads]);
        if r > maxr
            maxr = r;
            famInd = j;
            rawInd = ind;
        end
    end
    %     if maxr < LOW_ABUNDANCE_THRESHOLD
    %         raws = [bin(i).fam.raw];
    %         reals{i} = int2nt(consensusOut(scoreMat,{raws.seq},[raws.reads]));
    %     else
    reals{i} = int2nt(bin(i).fam(famInd).raw(rawInd).seq);
    %     end
end
reads = [bin.R]; %#ok<NASGU>
D = zeros(length(bin));
for b1 = 1:length(bin)
    for b2 = b1:length(bin)
        [~,al] = nwalign(reals{b1},reals{b2});
        D(b1,b2) = sum(al(1,:) ~= al(3,:));
        D(b2,b1) = D(b1,b2);
    end
end

err = zeros(4,4); %number of errors of each type
errCodon = zeros(4,4,4,4); %errors of each type in each context
bases = zeros(1,4); %total number of each base read
codon = zeros(4,4,4); %total number of each codon read

for b = 1:length(bin)
    bases = bases + bin(b).R * bin(b).baseCount;
    warning off Bioinfo:codoncount:UnknownSymbols
    warning off Bioinfo:codoncount:OtherSymbols
    [~,c1] = codoncount(reals{b},'frame',1);
    [~,c2] = codoncount(reals{b},'frame',2);
    [~,c3] = codoncount(reals{b},'frame',3);
    codon = codon + bin(b).R * (c1 + c2 + c3);
    t = bin(b).seq; %a template
    for f = 1:length(bin(b).fam)
        s = bin(b).fam(f).seq; %a sequence
        for i = find(t ~= s)
            if t(i) ~= 5 && s(i) ~= 5
                err(t(i),s(i)) = err(t(i),s(i)) + bin(b).fam(f).r;
                if i ~= 1 && i ~= length(t) && ...
                        t(i-1) ~= 5 && t(i+1) ~= 5
                    errCodon(t(i-1),t(i),s(i),t(i+1)) = ...
                        errCodon(t(i-1),t(i),s(i),t(i+1)) + ...
                        bin(b).fam(f).r;
                end
            end
        end
    end
end
save(fname,'bin','reads','reals','D','err','errCodon','bases','codon');
end

function lambda2 = lambdaShift(seq,subPos,subNT,lambda1,ERR,ERRa,context)
%computes lambda for a sequence to have been produced by a template based
%on the set of substitutions it has relative to that template
lambda2 = lambda1;
for i = 1:length(subPos)
    p = subPos(i);
    sub = subNT(i);
    if  context && p > 1 && p < length(seq)
        lambda2 = lambda2 * ERR(seq(p-1),seq(p),sub,seq(p+1)) / ...
            ERR(seq(p-1),seq(p),seq(p),seq(p+1));
    else
        lambda2 = lambda2 * ERRa(seq(p),sub) / ERRa(seq(p),seq(p));
    end
end
end

function [ERRa,ERR] = MLE(outFolder,iteration)
%compute the maximum likelihood error probabilities given the bin struct
%arrays output by a particular round of clustering
err = zeros(4,4); %number of errors of each type
errCodon = zeros(4,4,4,4); %errors of each type in each context
bases = zeros(1,4); %total number of each bases read
codon = zeros(4,4,4); %total number of each codon read

ERR = zeros(4,4,4,4);
ERRa = zeros(4,4);

files = dir([outFolder '/' num2str(iteration)]);
files = files([files.bytes]~=0); %remove empty els
for a = 1:length(files)
    s = load([outFolder '/' num2str(iteration) '/' files(a).name]);
    err = err + s.err;
    errCodon = errCodon + s.errCodon;
    bases = bases + s.bases;
    codon = codon + s.codon;
end

%add pseudocounts
err = err + 1;
errCodon = errCodon + 1;
bases = bases + 3; %num of context-indep errs added for each base
codon = codon + 12; %num of context-dep errs added for each context

for nt = 1:4
    err(nt,nt) = bases(nt) - sum(err(nt,setdiff(1:4,nt)));
    for l = 1:4
        for r = 1:4
            errCodon(l,nt,nt,r) = codon(l,nt,r) - ...
                sum(errCodon(l,nt,setdiff(1:4,nt),r));
        end
    end
end

for b = 1:4
    for e = 1:4
        ERRa(b,e) = err(b,e) / bases(b);
        for l = 1:4
            for r = 1:4
                ERR(l,b,e,r) = errCodon(l,b,e,r) / codon(l,b,r);
            end
        end
    end
end

save([outFolder '/' num2str(iteration) '/ERR'],'ERR','ERRa','err','bases');
end

function Cseq = consensus(scoreMat,seqs,reads)
%compute a consensus sequence for a particular cluster
%for now, just take the most abundance sequence
seqs = seqs(reads == max(reads));
reads = reads(reads == max(reads));
seqs = cellfun(@int2nt,seqs,'UniformOutput',false);
Cseq = seqs{1};
Cseq = uint8(nt2int(Cseq));
end

function Cseq = consensusOut(scoreMat,seqs,reads)
%compute a consensus sequence for a particular cluster
seqs = seqs(reads == max(reads));
reads = reads(reads == max(reads));
seqs = cellfun(@int2nt,seqs,'UniformOutput',false);
if length(seqs) == 1
    Cseq = seqs{1};
elseif length(seqs) == 2
    [~,al] = nwalign(seqs{1},seqs{2},'alphabet','nt',...
        'ScoringMatrix',scoreMat,'Glocal',true);
    Cseq = seqconsensus(al([1 3],:),...
        'ScoringMatrix',scoreMat,'alphabet','nt','gaps','all');
    Cseq( Cseq == '-' ) = [];
else
    m = cellstr(multialign(seqs,'ScoringMatrix',scoreMat,...
        'terminalGapAdjust',true));
    s = cell(1,sum(reads));
    for i = 1:length(reads);
        if i == 1
            s(1:reads(1)) = deal(m(i));
        else
            s(1+sum(reads(1:i-1)):sum(reads(1:i))) = deal(m(i));
        end
    end
    Cseq = seqconsensus(s,'ScoringMatrix',scoreMat,'alphabet','nt',...
        'gaps','all');
    Cseq(Cseq == '-') = [];
end
Cseq = uint8(nt2int(Cseq));
end