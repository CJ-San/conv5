//////////////////////////////////////////////////////////////////////
`timescale 1ps/1ps
module tb;
//////////////////////////parmaeters////////////////////////////////////
localparam	CLK_PERIOD=1031;
parameter	RUN_TIME=67051085;
//////////////////////////////tb_interface/////////////////////////////
logic	[127:0]	image_r	[65024:0];
logic	[127:0]	image_g	[65024:0];
logic	[127:0]	image_b	[65024:0];
logic	[7:0]	results_mem [65024:0];
logic	[7:0]	golden_results [65024:0];
logic	signed [71:0]	filter_r;
logic	signed [71:0]	filter_g;
logic	signed [71:0]	filter_b;
int i,j;

///////////////////////////DUT interface////////////////////////////////////////////////////
logic	[127:0]	dut_image_r;
logic	[127:0]	dut_image_g;
logic	[127:0]	dut_image_b;
logic	[(8*9)-1:0]	dut_kernel_r;
logic	[(8*9)-1:0]	dut_kernel_g;
logic	[(8*9)-1:0]	dut_kernel_b;
logic	[15:0]	dut_im_address_r;
logic	[15:0]	dut_im_address_g;
logic	[15:0]	dut_im_address_b;
logic	[15:0]	dut_res_address;
logic	[7:0]	dut_results;
logic 	dut_give_im_r,dut_give_im_g,dut_give_im_b,clk_tb,rst,dut_take_results;
//////////////////////////////////Instance of DUT///////////////////////////////////////
conv_pool dut(	.clk(clk_tb),
				.rst(rst),
				.image_4x4_r(dut_image_r),
				.image_4x4_g(dut_image_g),
				.image_4x4_b(dut_image_b),
				.kernel_r(dut_kernel_r),
				.kernel_g(dut_kernel_g),
				.kernel_b(dut_kernel_b),
				.input_re(dut_give_im_r),
				.input_addr(dut_im_address_r),
				.output_we(dut_take_results),
				.output_addr(dut_res_address),
				.y(dut_results)
				);
assign dut_give_im_g=dut_give_im_r;
assign dut_give_im_b=dut_give_im_r;
assign dut_im_address_g=dut_im_address_r;
assign dut_im_address_b=dut_im_address_r;
///////////////////////////clk and rst////////////////////////////////////////////////////
always 
	begin 
		#(CLK_PERIOD/2) clk_tb=~clk_tb; 
	end
task reset_dut();
  #(CLK_PERIOD/2) rst=0;
  @(posedge clk_tb);
     #(CLK_PERIOD/3)  rst=1;	
endtask
/////////////////////////initalize image,filters////////////////////////////////////
task initialize_image();
		integer fp1,fp2,fp3,p;
		fp1=$fopen("./image_r.txt","r");
		for(i=0;i<65025;i++)
			  p=$fscanf(fp1,"%h\n",image_r[i]);
		$fclose(fp1);        
		fp2=$fopen("./image_g.txt","r");
		for(i=0;i<65025;i++)
			  p=$fscanf(fp2,"%h\n",image_g[i]);
		$fclose(fp2);        
		fp3=$fopen("./image_b.txt","r");
		for(i=0;i<65025;i++)
			  p=$fscanf(fp3,"%h\n",image_b[i]);
		$fclose(fp3);
		        
endtask
task initialize_filter();
		integer fp1,fp2,fp3,p;
		fp1=$fopen("./filter_r.txt","r");
			  p=$fscanf(fp1,"%h\n",filter_r);
		$fclose(fp1);        
		fp2=$fopen("./filter_g.txt","r");
			  p=$fscanf(fp2,"%h\n",filter_g);
		$fclose(fp2);        
		fp3=$fopen("./filter_b.txt","r");
			  p=$fscanf(fp3,"%h\n",filter_b);
		$fclose(fp3);        
endtask
task initialize_gloden_results();
integer fp1,p;
  fp1=$fopen("./golden_results.txt","r");
		for(i=0;i<65025;i++)
			p=$fscanf(fp1,"%h\n",golden_results[i]);
		$fclose(fp1);  
endtask

////////////////////////////////filters/////////////////////////////////////////////////////////////////
assign dut_kernel_r={filter_r};
assign dut_kernel_g={filter_g};
assign dut_kernel_b={filter_b};
////////////////////////////memory_read_models///////////////////////////////////////////////////////////
always @(posedge clk_tb)
begin
	if	(dut_give_im_r==1) dut_image_r<=image_r[dut_im_address_r];
	else dut_image_r<=0;
end
always @(posedge clk_tb)
begin
	if	(dut_give_im_g==1) dut_image_g<=image_g[dut_im_address_g];
	else dut_image_g<=0;
end
always @(posedge clk_tb)
begin
	 if	(dut_give_im_b==1) dut_image_b<=image_b[dut_im_address_b];
	else dut_image_b<=0;
end
////////////////////////////memory_write_models///////////////////////////////////////////////////////////
//initial results_mem=0;
always @(posedge clk_tb)
begin
	
	if	(dut_take_results==1) results_mem[dut_res_address]<=dut_results;
end
//////////////////////////////////DUT verifier////////////////////////////////////////////////////////////////////
task verify_dut();
	bit error_check;
	int i;
	error_check=0;
	for(i=0;i<65025;i++)
		begin
			
			if((results_mem[i])!==(golden_results[i]))
			begin
				error_check=1;
				$display("FAILED...!!! Check results of block %d,result:%d,golden:%d\n",i,results_mem[i],golden_results[i]);
				break;
			end
		end
	if(error_check==0)$display("PASS..!!\n");
endtask
//////////////////////////////////////////////////////////////////////////////////////////////////////
initial 
	begin
		
      	clk_tb=0;
		initialize_filter();
		initialize_gloden_results();
		initialize_image();
		reset_dut();
      		#(RUN_TIME);
		verify_dut();
		$finish();
	end	
/*initial 
begin   
	 $dumpfile("./convpool_output.vcd");
	 $dumpvars (0,tb);
end*/
/*initial begin
	$shm_open("waves.shm");
	$shm_probe(tb, "ascm");
end*/
endmodule	
