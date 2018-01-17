function plotMotionFromTxt(allPath);
% Convert all '*rp_txt' to motion plots(.png) in the given path.
    filt = ['^rp_*','.*\.txt$'];
    %    b = spm_select([Inf],'any','Select realignment parameters',[],pwd,filt);
       [files,dirs] = spm_select('FPListRec',allPath,filt);
        b=files;

       for i = 1:size(b,1)

        [p nm e v] = spm_fileparts(b(i,:));

        printfig = figure;
        set(printfig, 'Name', ['Motion parameters: subject ' num2str(i) ],'Visible', 'on');
        loadmot = load(deblank(b(i,:)));
        subplot(2,1,1);
        plot(loadmot(:,1:3));
        grid on;
        % ylim(scaleme);  % enable to always scale between fixed values
        title(['Motion parameters: shifts (top, in mm) and rotations (bottom,in dg)'], 'interpreter', 'none');

        subplot(2,1,2);
        plot(loadmot(:,4:6)*180/pi);
        grid on;
        %  ylim(scaleme);
        title(['Data from ' p], 'interpreter', 'none');
        filename = char(regexprep((regexp(b(i,:),'(?=df\d+).*(?=.txt)','match')),'\','_'));

        path=char(regexp(b(i,:),'.*(?=\\rp)','match'));
        motname = [path filesep 'motion_sub_' filename '.png'];
        print(printfig, '-dpng', '-noui', '-r100', motname);
        close(printfig)
       end
end