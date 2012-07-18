function test_suite = testComplex
initTestSuite;

    % BEGIN INITIALIZATION

    function opt = setup
        % Make sure environment is pristine
        clear all
        clear functions
        
        % Make sure SPGL1 exists in path
        addpath('../..')
        
        % Fix random generator for repeatable experiments
        defaultStream = RandStream.getDefaultStream;
        savedState = defaultStream.State;
        RandStream.setDefaultStream(RandStream('mt19937ar','seed',8888));
        
        % Make test cases
        m = 3000;
        n = 10000;
        k = 300;
        
        A = randn(m,n) + i*randn(m,n);
        x = zeros(n,1);
        x(floor((n-1)*rand(k,1))+1) = randn(k,1) + i*randn(k,1); % inject k-sparse signal in random places
        b = A*x;
        
        % Set sparse recovery element-wise relative tolerance
        recov_tol = 2e-5;
        
        % Inject options into test cases
        opt.A = A;
        opt.b = b;
        opt.x = x;
        opt.recov_tol = recov_tol;
        opt.savedState = savedState;

    function teardown(opt)
        % Restore random stream and restore variables
        defaultStream.State = opt.savedState;
        clear opt
        
        % Remove the test SPGL1 from path
        rmpath('../..')
        clear all
        clear functions
    


    % BEGIN TESTS
    
    function testStandardCplxSPGL1(opt)
        % figure;plot(opt.x)
        spg_opts  = spgSetParms('verbosity' ,1);
        [x, r, g, info] = spgl1(opt.A, opt.b, [], 1e-2, [], spg_opts);
        % figure;plot(x)
        % figure;plot(x-opt.x)
        norm(x-opt.x,inf)
        assertElementsAlmostEqual(x, opt.x, 'absolute', opt.recov_tol);
        assertTrue(info.iter < 71);
