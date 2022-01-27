function t2 = interface()
    disp('hello python')
    [i,j, k] = foobar();
    t(1).foo = i;
    t(1).bar = j; 
    t(1).num = k;

    [i,j, k] = foobar2();
    t(2).foo = i;
    t(2).bar = j; 
    t(2).num = k;
    
    t2.one = t(1);
    t2.two{1} = t(2);
    t2.two{2} = t(1);
end

function [foo, bar, num] = foobar()
    foo = 'foo';
    bar = 'bar';
    num = 3;
end

function [foo, bar, num] = foobar2()
    foo = 'foo2';
    bar = 'bar2';
    num = 32;
end