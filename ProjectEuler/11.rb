#Recursive function that will return the product of n consecutive elements in the grid
def grid_product(row=0, col=0, n=4, dir="")
  if n==0 || row >= $nRows || col >= $nCols then
    return 1
  else
    case dir
    when "south"
      return $grid[row][col] * grid_product(row+1, col, n-1, dir)
    when "east"
      return $grid[row][col] * grid_product(row, col+1, n-1, dir)
    when "northeast"
      return $grid[row][col] * grid_product(row-1, col+1, n-1, dir)
    when "southeast"
      return $grid[row][col] * grid_product(row+1, col+1, n-1, dir)
    when "southwest"
      return $grid[row][col] * grid_product(row+1, col-1, n-1, dir)
    end
  end
end

#Open and translate file contents to a 2D array
$grid = []
File.open('11.txt','r') do |file|
  file.each_line do |line|
    #line.chomp!
    nums = line.chomp.split(' ').collect {|num_str| num_str.to_i}
    $grid << nums
  end
end

#Iterate through and find max possible products at each index
$max = 0
$dir = ["south", "east", "northeast", "southeast", "southwest"]
$nRows = $grid.size
$nCols = $grid[0].size

(0...$nRows).each do |i|
  (0...$nCols).each do |j|
    results = []
    $dir.each do |direction|
      results << grid_product(i,j,4,direction)
    end

    if results.max > $max then
      $max = results.max
    end
  end
end

puts "Answer is #{$max}"
