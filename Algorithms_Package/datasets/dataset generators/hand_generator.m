function [] = hand_generator(digit1,digit2,digit3,comb)
    % Digit1,2,3 are the hand written digits that will
    % be used in the dataset that will be generated. 
    % They start from 0 to 9. The whole input parameter
    % takes 2, 1 or 0 as values and defines if we will use 
    % both training and test set (2), just the training set (1)
    % or just the test set (0) of the hand_written_digits8x8
    % dataset.
    
    load hand_written_digits8x8;
    cd('datasets/');
    str0 = 'hand_written_digits8x8_3of10_';
    str1 = num2str(digit1);
    str2 = num2str(digit2);
    str3 = num2str(digit3);
    str4 = 'small.mat';
    str5 = 'medium.mat';
    str6 = 'large.mat';
    
    switch comb
        case 0
            Ct1 = Ctest(digit1*200+1:digit1*200+200,:);
            Ct2 = Ctest(digit2*200+1:digit2*200+200,:);
            Ct3 = Ctest(digit3*200+1:digit3*200+200,:);
            C = [Ct1;Ct2;Ct3];
            Xt1 = Xtest(digit1*200+1:digit1*200+200,:);
            Xt2 = Xtest(digit2*200+1:digit2*200+200,:);
            Xt3 = Xtest(digit3*200+1:digit3*200+200,:);
            X = [Xt1;Xt2;Xt3];
            strOut = strcat(str0,str1,str2,str3,str4);
            save(strOut,'X','C'); 
            fprintf('The %s dataset was created\n\n',strOut);
        case 1
            Ctr1 = Ctrain(digit1*700+1:digit1*700+700,:);
            Ctr2 = Ctrain(digit2*700+1:digit2*700+700,:);
            Ctr3 = Ctrain(digit3*700+1:digit3*700+700,:);
            C = [Ctr1;Ctr2;Ctr3];
            Xtr1 = Xtrain(digit1*700+1:digit1*700+700,:);
            Xtr2 = Xtrain(digit2*700+1:digit2*700+700,:);
            Xtr3 = Xtrain(digit3*700+1:digit3*700+700,:);
            X=[Xtr1;Xtr2;Xtr3];
            strOut = strcat(str0,str1,str2,str3,str5);
            save(strOut,'X','C');
            fprintf('The %s dataset was created\n\n',strOut);
        case 2
            Ct1 = Ctest(digit1*200+1:digit1*200+200,:);
            Ct2 = Ctest(digit2*200+1:digit2*200+200,:);
            Ct3 = Ctest(digit3*200+1:digit3*200+200,:);
            Ctr1 = Ctrain(digit1*700+1:digit1*700+700,:);
            Ctr2 = Ctrain(digit2*700+1:digit2*700+700,:);
            Ctr3 = Ctrain(digit3*700+1:digit3*700+700,:);
            C = [Ct1;Ctr1;Ct2;Ctr2;Ct3;Ctr3];
            Xt1 = Xtest(digit1*200+1:digit1*200+200,:);
            Xt2 = Xtest(digit2*200+1:digit2*200+200,:);
            Xt3 = Xtest(digit3*200+1:digit3*200+200,:);
            Xtr1 = Xtrain(digit1*700+1:digit1*700+700,:);
            Xtr2 = Xtrain(digit2*700+1:digit2*700+700,:);
            Xtr3 = Xtrain(digit3*700+1:digit3*700+700,:);
            X = [Xt1;Xtr1;Xt2;Xtr2;Xt3;Xtr3];
            strOut = strcat(str0,str1,str2,str3,str6);
            save(strOut,'X','C');
            fprintf('The %s dataset was created\n\n',strOut);
        otherwise
            fprintf('comb needs a value of 0, 1 or 2 to work\n\n');
    end
    cd('..');
end